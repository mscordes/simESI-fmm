#include "State.h"
#include "Atom.h"
#include "Constants.h"
#include "StringUtils.h"

#include <iomanip>
#include <fstream>
#include <stdexcept>

// Implementation for building the state from a protein .pdb file
// NOTE: This is only intended to parse protein specific atoms
void State::parseProteinPDB(const string& pdbFile) {
    proteinMap.clear();
    residueSet.clear();

    ifstream file(pdbFile);
    if (!file.is_open()) throw runtime_error("Error opening .pdb file: " + pdbFile);

    int resCount = -1;
    string currentResID = "";
    shared_ptr<Residue> residue;

	int chainCount = -1;
    string currentChain = "";
    char chain{};
    shared_ptr<Protein> protein;

    string line;
    while (getline(file, line)) {

        if (line.substr(0, 4) != "ATOM") continue;

        // Get residue ID and chain
		string resName = trim(line.substr(17, 3)); // NOTE: This trims resName to only 3 characters!
        int tempResNum = stoi(line.substr(22, 4)); // resNum from .pdb not always equal to resCount
		string tempResID = to_string(tempResNum) + resName;

        // Cycle residues
        if (tempResID != currentResID) {

            // Check if amino acid residue
            if (Constants::aminoAcids.find(resName) == Constants::aminoAcids.end()) continue;

            // Store previous residue
            if (currentResID != "") {
                residueSet[resName].insert(residue);
                protein->addResidue(residue);
            }
            currentResID = tempResID;

            // Cycle chains
            string tempChain = trim(line.substr(21, 1)); // Chain ID from .pdb
            if (tempChain != currentChain) {

                if (currentChain != "") proteinMap[chain] = protein;
                currentChain = tempChain;

                ++chainCount;
                if (chainCount > 25) throw runtime_error("Chain limit exceeded 25 during parseProteinPDB().");

                chain = Constants::alphabet.at(chainCount);
                protein = make_shared<Protein>(chain);
			}

			// Create new residue
            ++resCount;
            residue = make_shared<Residue>(
				to_string(resCount) + resName, // Correct residue ID
				resName,                       // Residue name
				resCount,                      // Correct residue number
				chain                          // Chain
            );
		}

        // Build atom
        string atomName = trim(line.substr(12, 4));
		residue->addAtom(
            make_shared<Atom>(
                atomName,               // Atom name 
                atomName[0],            // Element
                Vec3D(                  // Coordinates (in nm)
                    stod(line.substr(30, 8)) / 10,
                    stod(line.substr(38, 8)) / 10,
                    stod(line.substr(46, 8)) / 10),
                Vec3D(0.0, 0.0, 0.0),   // Velocity (pdbs dont contain vel info)
                residue                 // Parent residue
            )
        );
        
    }
    file.close();

    residueSet[residue->name].insert(residue);
    protein->addResidue(residue);
    proteinMap[chain] = protein;

    // Set GPB's
    for (auto& [resType, residues] : residueSet) {
        if (Constants::GPBs.contains(resType)) {
            double gpb = Constants::GPBs.at(resType);
            for (auto& residue : residues) residue->gpb = gpb;
        }
    }

    // Set termini info 
    for (auto& [chain, protein] : proteinMap) {
        auto& nterm = protein->residues.front();
        nterm->termini = true;
        nterm->gpb = Constants::GPBs.at("NTERM");
        auto& cterm = protein->residues.back();
        cterm->termini = true;
        cterm->gpb = Constants::GPBs.at("CTERM");
    }
}