#include "TopologyUtils.h"
#include "Atom.h"
#include "Constants.h"
#include "FileUtils.h"
#include "Parameters.h"
#include "Protein.h"
#include "Residue.h"
#include "RNG.h"
#include "Subprocess.h"
#include "StringUtils.h"

#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <stdexcept>
#include <sstream>

// Termini codes for pdb2gmx change depending on amino acid type, account for that here
// Return pair of codes where first = prot, second = deprot
pair<string, string> getNtermPdb2GmxCodes(const shared_ptr<Residue>& residue) {
    string resType = residue->name;
    if (resType == "MET") return { "1", "2" };
    else if (resType == "PRO") return { "1", "0" };
    else return { "0", "1" };
}

// Helper function for setPdb2GmxInputs()
// Updates in place pdb2gmx inputs with the pdb2gmx code for residue given an inputted pH and residue pKa
static void setProtState(
    const shared_ptr<Residue>& residue, 
    double pH, 
    vector<string>& chainInputs, 
    bool isNterm, 
    bool isCterm) 
{
    double pka;
    pair<string, string> codes;

    try {
        if (isNterm) {
            pka = residue->termini_pka.value();
            codes = getNtermPdb2GmxCodes(residue);
        }
        else if (isCterm) {
            pka = residue->termini_pka.value();
            codes = Constants::pdb2gmxCodes.at("CTERM");    
        }
        else {
            pka = residue->pka.value();                 
            codes = Constants::pdb2gmxCodes.at(residue->name); 
        }
    }
    catch (const bad_optional_access& e) {
        throw runtime_error(
            "Could not find pKa value for residue " + residue->ID +
            " (missing optional): " + e.what());
    }
    catch (const out_of_range& e) {
        throw runtime_error(
            "Could not find pdb2gmx code for residue " + residue->ID +
            " (map lookup failed): " + e.what());
    }

    double prob = 1.0 / (1 + pow(10.0, pH - pka)); // Henderson Hasselbach
    bernoulli_distribution d(prob);
    const bool isDeprotonated = bernoulli_distribution(prob)(RNG::gen);
    chainInputs.push_back(isDeprotonated ? codes.first : codes.second);
}

// If an amino acid is protonated
bool isAA_Protonated(const shared_ptr<Residue>& residue, bool isNterm, bool isCterm) {
    string resName;
    if (isNterm) resName = "NTERM";
    else if (isCterm) resName = "CTERM";
    else resName = residue->name;

    int numHydrogens = 0;
    const auto& hydrogens = Constants::hydrogensMap.at(resName);

    for (const auto& atom : residue->atoms) 
        if (hydrogens.contains(atom->name)) 
            ++numHydrogens;
    return (hydrogens.size() == numHydrogens); 
}

// If deprotonated histidine HISD, else return false for HISE
bool isHISD(const shared_ptr<Residue>& his) {
    for (const auto& atom : his->atoms) 
        if (atom->name == "HD1") return true;
    return false;
}

// Set inputs for pdb2gmx probabilistically relative to a given pH
map<char, vector<string>> setPdb2GmxInputs(State& state, double ph) {
    const auto& proteinMap = state.proteinMap;
    const auto& residueSet = state.residueSet;

    map<char, vector<string>> inputs;
    for (const auto& [chain, protein] : proteinMap) {
        vector<string> chainInputs;

        for (const auto& resType : Constants::pdb2gmxOrder) 
            for (const auto& residue : protein->residues) 
                if (residue->name == resType) 
                    setProtState(residue, ph, chainInputs, false, false);

        const auto& nterm = protein->residues.front();
        setProtState(nterm, ph, chainInputs, true, false);

        const auto& cterm = protein->residues.back();
        setProtState(cterm, ph, chainInputs, false, true);

        inputs[chain] = chainInputs;
    }
    return inputs;
}

// Get protonation state and corresponding pdb2gmx input of a residue
static void getProtState(
    const shared_ptr<Residue>& residue, 
    vector<string>& chainInputs, 
    bool isNterm, 
    bool isCterm) 
{
    pair<string, string> codes;

    if (isNterm) {
        codes = getNtermPdb2GmxCodes(residue);
        chainInputs.push_back(isAA_Protonated(residue, true, false) ? codes.first : codes.second);
    }
    else if (isCterm) {
        auto it = Constants::pdb2gmxCodes.find("CTERM");
        if (it != Constants::pdb2gmxCodes.end()) 
            codes = it->second;
        else 
            throw runtime_error("Missing CTERM entry in pdb2gmxCodes.");

        chainInputs.push_back(isAA_Protonated(residue, false, true) ? codes.first : codes.second);
    }
    else if (residue->name == "HIS") {
        if (isAA_Protonated(residue))
            chainInputs.push_back("2");
        else if (isHISD(residue)) 
            chainInputs.push_back("0");
        else
            chainInputs.push_back("1");
    }
    else {
        auto it = Constants::pdb2gmxCodes.find(residue->name);
        if (it != Constants::pdb2gmxCodes.end()) 
            codes = it->second;
        else
            throw runtime_error("Missing entry for residue: " + residue->name);
        chainInputs.push_back(isAA_Protonated(residue) ? codes.first : codes.second);
    }
}

// Returns pdb2gmx inputs to maintain protonation states of State
map<char, vector<string>> getPdb2GmxInputs(State& state) {
    const auto& proteinMap = state.proteinMap;
    const auto& residueSet = state.residueSet;

    map<char, vector<string>> inputs;
    for (const auto& [chain, protein] : proteinMap) {
        vector<string> chainInputs;

        for (const auto& resType : Constants::pdb2gmxOrder)
            for (const auto& residue : protein->residues)
                if (residue->name == resType)
                    getProtState(residue, chainInputs, false, false);

        const auto& nterm = protein->residues.front();
        getProtState(nterm, chainInputs, true, false);

        const auto& cterm = protein->residues.back();
        getProtState(cterm, chainInputs, false, true);

        inputs[chain] = chainInputs;
    }
    return inputs;
}

// Convert .top file of a given protein chain to an .itp file
static void convertTop2Itp(const string& topFile, char chain) {
    if (!fs::exists(topFile))
        throw runtime_error("Could not find the .top file " + topFile + " for protein chain"
            ".\nLikely caused by failure of pdb2gmx due to improper input structure.");

    ifstream inFile(topFile);
    if (!inFile.is_open())
        throw runtime_error("Could not open the .top file " + topFile + ".");

    string itpFile = fs::path(topFile).stem().string() + ".itp";
    ofstream outFile(itpFile);
    if (!outFile.is_open()) 
        throw runtime_error("Could not create .itp file " + itpFile + ".");

    string line;
    bool writeLine = false;
    while (getline(inFile, line)) {
        if (line.find("[ moleculetype ]") != string::npos) writeLine = true;
        if (line.find("; Include Position restraint file") != string::npos) break;

        if (line.find("Protein") != string::npos) line = string(1, chain) + "          3";
        if (writeLine) outFile << line << "\n";
    }
    inFile.close();
    outFile.close();
}

// Make .top file with proper .itp info and residue ordering/numbering
void createTOP(const string& fname, const State& state) {
    ofstream file(fname);
    if (!file.is_open()) throw runtime_error("Could not create .top file " + fname + ".");

    vector<string> chains;
    for (const auto& [chain, protein] : state.proteinMap) {
        chains.push_back(string(1, chain));
    }

    file << "#include \"charmm36.ff/forcefield.itp\"\n";
    for (const auto& chain : chains) file << "#include \"" << chain << ".itp\"\n";
    file << "#include \"o2.itp\"\n";         
    file << "#include \"n2.itp\"\n";           
    file << "#include \"nh4.itp\"\n";         
    file << "#include \"nh3.itp\"\n";          
    file << "#include \"aceh.itp\"\n";       
    file << "#include \"ace.itp\"\n";         
    file << "#include \"oh.itp\"\n";    
    file << "#include \"h3o.itp\"\n";    
    file << "#include \"tip4p_2005.itp\"\n\n";  
    file << "[ system ]\n";
    file << "simESI System\n\n";
    file << "[ molecules ]\n";
    for (const auto& chain : chains) file << chain << "                   " << "1\n";

    for (const auto& resType : Constants::topOrder) {
        if (auto it = state.residueSet.find(resType); it != state.residueSet.end())
            file << resType << "                 " << it->second.size() << '\n';
    }

    file.close();
}

// Call pdb2gmx and create .gro, .top and .itp file(s) from given protonation states
void call_pdb2gmx(const string& topName, const map<char, vector<string>>& gmxInputs, State& state, bool keepGro) {

    static const double boxD = Parameters::Get().getBoxSize();
    static const string& gmxEnv = Parameters::Get().getGMXEnv();

    // Create unique .itp/.gro files for each monomer
    for (const auto& [chain, protein] : state.proteinMap) {
        string chainStr = string(1, chain);

        string newGroName = "top_" + chainStr + ".gro";
        string preGroName = "pre_" + chainStr + ".gro";

        deleteFile(preGroName);
        deleteFile(newGroName);
        state.writeChainGRO(preGroName, boxD, *protein);

        { // Make chain specific .top file
            if (gmxInputs.find(chain) == gmxInputs.end()) 
                throw runtime_error("Chain " + chainStr + " not found in pdb2gmx inputs.");
            const vector<string>& chainInputs = gmxInputs.at(chain);

            string newTopName = chainStr + ".top";
            deleteFile(newTopName);
    

            string command = gmxEnv + " pdb2gmx -f " + preGroName + " -o " + newGroName + " -p " + newTopName
                + " -ff charmm36 -water none -lys -arg -glu -asp -his -ter -ignh";
            subprocess(command, chainInputs);

            if (!fs::exists(newTopName))
                throw runtime_error(".top file not created for protein chain " + chainStr + ".");

            // Convert .top to .itp
            convertTop2Itp(newTopName, chain);

            deleteFile(fs::path(newTopName));
            deleteFile(fs::path(preGroName));
            deleteFile(fs::path("posre.itp"));
        }
    }

    // If keeping .gro need to stitch chains together
    if (keepGro) {

        string groFile = fs::path(topName).stem().string() + ".gro";
        ofstream outFile(groFile);
        if (!outFile.is_open()) throw runtime_error("Could not create .gro file " + groFile + ".");

        outFile << "simESI-fmm generated .gro file\n";

        // We don't know atom count yet, so place hold for now
        streampos countPos = outFile.tellp();
        outFile << setw(10) << "\n";

        int atomCount = 0;
        for (const auto& [chain, protein] : state.proteinMap) {
            string chainGroFile = "top_" + string(1, chain) + ".gro";
            if (!fs::exists(chainGroFile))
                throw runtime_error("Could not find the .gro file " + chainGroFile + " for protein chain.");

            ifstream inFile(chainGroFile);
            if (!inFile.is_open())
                throw runtime_error("Could not open the chain .gro file " + chainGroFile + ".");

            string line; // Skip headers
            getline(inFile, line);
            getline(inFile, line);

            while (getline(inFile, line)) {

                // Ignore last line with box dimensions
                if (line.size() < 44) break;

                // Correct atom count
                ++atomCount;
                ostringstream oss; 
                oss << setw(5) << right << atomCount; 
                string countField = oss.str();
                line.replace(15, 5, countField);

                outFile << line << "\n";
            }
            inFile.close();
        }

        outFile << fixed << setprecision(5) << "  " << boxD << "  " << boxD << "  " << boxD << '\n';

        // Edit number of atoms in second line
        outFile.seekp(countPos);
        outFile << atomCount;

        outFile.close();
    }

    // Clean up unneeded chain .gro files
    for (const auto& [chain, protein] : state.proteinMap)
        deleteFile(fs::path(string(1, chain) + ".gro"));

    // Create the master .top file
    createTOP(topName, state);
}

static vector<pair<string, string>> parseITP(const string& fname) {
    ifstream f(fname);
    if (!f.is_open()) 
        throw runtime_error("Error opening .itp file: " + fname);

    vector<pair<string, string>> atoms;
    string line;
    while (getline(f, line)) {
        vector<string> words = splitLine(line);
        if (!line.empty() && line.find_first_not_of(' ') != string::npos) {
            if (words[1] == "bonds") break;
            if (words[0][0] != ';' && words[0][0] != '[' && words.size() >= 7) 
                atoms.emplace_back(words[3].substr(0, 3), words[4]);
        }
    }
    return atoms;
}

static unordered_map<string, int> parseTOP(const string& fname) {
    ifstream f(fname);
    if (!f.is_open())
        throw runtime_error("Error opening .top file: " + fname);
    const auto& topOrder = Constants::topOrder; 

    unordered_map<string, int> resMap;
    string line;
    while (getline(f, line)) {
        vector<string> words = splitLine(line);
        if (words.size() == 2 &&
            find(topOrder.begin(), topOrder.end(), words[0]) != Constants::topOrder.end())
        {
            resMap[words[0]] = stoi(words[1]);
        }
    }
    return resMap;
}

// Test for congruency between State and .top
void testTOP(const State& s, string top) {
    auto itpResMap = parseTOP(top);
    for (const auto& [resType, residues] : s.residueSet) {
        if (Constants::aminoAcids.contains(resType)) continue;

        auto it = itpResMap.find(resType);
        if (it == itpResMap.end())
            throw runtime_error("Could not find State resType " + resType + " in top file.");
        auto itpCount = it->second;

        int stateCount = static_cast<int>(residues.size());
        if (itpCount != stateCount) 
            throw runtime_error("Residue count mismatch for " + resType + ". " + 
                to_string(stateCount) + " in state versus " + to_string(itpCount) + " in top.");
    }
    cout << "Matched residue counts in State and top.\n";

    vector<shared_ptr<Atom>> stateProtAtoms;
    for (const auto& [chain, p] : s.proteinMap)
        for (const auto& r : p->residues)
            for (const auto& a : r->atoms)
                stateProtAtoms.push_back(a);

    vector<pair<string, string>> itpProtAtoms;
    for (const auto& pair : s.proteinMap) {
        auto cAtoms = parseITP(string(1, pair.first) + ".itp");
        itpProtAtoms.insert(itpProtAtoms.end(), cAtoms.begin(), cAtoms.end());
    }

    size_t stateSize = stateProtAtoms.size();
    size_t itpSize = itpProtAtoms.size();
    if (stateSize > itpSize) {
        cout << "More atoms in State (" + to_string(stateSize) + ") than top (" 
            + to_string(itpSize) + ")\n";
    }
    else if (stateSize < itpSize) {
        cout << "Less atoms in State (" + to_string(stateSize) + ") than top ("
            + to_string(itpSize) + ")\n";
    }
    else {
        cout << "Atom count in State and .itp equal.";
    }
    size_t size = min(stateSize, itpSize);

    for (size_t i = 0; i < size; ++i) {
        const auto& sa = stateProtAtoms[i];
        const auto& ia = itpProtAtoms[i];

        if (sa->parent.lock()->name != ia.first) {
            string e = "Mismatch residue type of State (" + sa->parent.lock()->name +
                ") versus top (" + ia.first + ")\n";
            throw runtime_error(e);
        }

        if (sa->name != ia.second) {
            string e = "Mismatch atom at index of " + to_string(i) + " type of State (" + sa->name +
                " " + sa->parent.lock()->name + ") versus top (" + ia.second + " " + ia.first + ")\n";
            throw runtime_error(e);
        }   
    }
}