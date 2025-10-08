#include "ParsePKA.h"
#include "Atom.h"
#include "Constants.h"
#include "StringUtils.h"
#include "Subprocess.h"

#include <fstream>

// Update pkaMap of chain : (resID : pKa value) for all titratable residues in .pka file
void parsePKA(State& state) {
    unordered_map<char, unordered_map<string, double>> pkaMap;

    // Temporarily rename C-terminal OT2 to OXT for PROPKA
    for (const auto& [chain, protein] : state.proteinMap) 
        for (auto& atom : protein->residues.back()->atoms) 
            if (atom->name == "OT2") atom->name = "OXT";

    state.writeProteinPDB("pka.pdb");
    for (const auto& [chain, protein] : state.proteinMap) 
        for (auto& atom : protein->residues.back()->atoms) 
            if (atom->name == "OXT") atom->name = "OT2";

    // Call PROPKA3
    string command = "python -m propka pka.pdb";
    subprocess(command);

    // Parse .pka file into map of  chain : (resID : pKa value) for all titratable residues in .pka file
    ifstream inputFile("pka.pka");
    if (!inputFile.is_open())
        throw runtime_error("pKa file not found. Likely do not have PROPKA3 installed and callable from command line.");

    string line;
    while (getline(inputFile, line))
    {
        vector<string> words = splitLine(line);

        // Expect something like: <ResidueName> <ResidueNumber> <Chain> ... pKa ...
        if (words.size() == 5 && isInteger(words[1])) {

            char chain = words[2][0];
            int resNum = stoi(words[1]) - 1; // Convert 0-based

            // Update resNum to account for total res count
            try {
                const auto& chainResidues = state.proteinMap.at(chain)->residues;
                auto& sp = chainResidues.at(resNum);
                resNum = sp->num;
            }
            catch (...) {
                throw runtime_error("Out of bounds chain or residue when assigning pKa values.");
            }

            string resName = words[0];
            string resID = to_string(resNum) + resName;
            double pka = stod(words[3]);

            // Termini special naming handling
            if (resName == "N+" || resName == "C-") {
                resID = resName;
            }

            // Insert into nested map
            pkaMap[chain][resID] = pka;
        }
    }
    inputFile.close();

    // Assign pKa values to residues and their termini
    for (auto& [chain, protein] : state.proteinMap) {
        const auto& chainPkaMap = pkaMap.at(chain); 

        for (auto& r: protein->residues) {
            if (Constants::titAminoAcids.contains(r->name)) { 
                try { 
                    r->pka = chainPkaMap.at(r->ID); 
                } 
                catch (const out_of_range&) { 
                    if (Constants::defaultPkaVals.contains(r->name)) { 
                        r->pka = Constants::defaultPkaVals.at(r->name); 
                        cout << "WARNING: Could not find a PROPKA pKa value for residue " << r->ID << 
                            ". Using default pKa value of " << *r->pka << " instead.\n"; 
                    } 
                    else { 
                        throw runtime_error("Could not assign pKa value to residue " + r->ID + ".\n"); 
                    } 
                } 
            } 
        } 

        // Termini 
        auto& nterm = protein->residues.front(); 
        try { 
            nterm->termini_pka = chainPkaMap.at("N+"); 
        } 
        catch (const out_of_range&) { 
            if (Constants::defaultPkaVals.contains("N+")) { 
                nterm->termini_pka = Constants::defaultPkaVals.at("N+"); 
                cout << "WARNING: Could not find a PROPKA pKa value for N-terminal residue " << nterm->ID << 
                    ". Using default pKa value of " << *nterm->termini_pka << " instead.\n"; 
            } 
            else { 
                throw runtime_error("Could not assign pKa value to N-terminal residue " + nterm->ID + ".\n"); 
            } 
        } 

        auto& cterm = protein->residues.back(); 
        try { 
            cterm->termini_pka = chainPkaMap.at("C-"); 
        } 
        catch (const out_of_range&) { 
            if (Constants::defaultPkaVals.contains("C-")) { 
                cterm->termini_pka = Constants::defaultPkaVals.at("C-"); 
                cout << "WARNING: Could not find a PROPKA pKa value for C-terminal residue " << cterm->ID << 
                    ". Using default pKa value of " << *cterm->termini_pka << " instead.\n"; 
            } 
            else { 
                throw runtime_error("Could not assign pKa value to C-terminal residue " + cterm->ID + ".\n"); 
            } 
        } 
    }

    // Set other non-protein pKas
    for (const auto& resType : Constants::nonProtTitResidues) {
        if (state.residueSet.contains(resType)) {
            double pka = Constants::defaultPkaVals.at(resType);
            for (auto& residue : state.residueSet.at(resType))
                residue->pka = pka;
        }
    }
}