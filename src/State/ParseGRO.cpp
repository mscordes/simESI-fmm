#include "State.h"
#include "Atom.h"
#include "Constants.h"
#include "CoordinateManip.h"
#include "StringUtils.h"

#include <iomanip>
#include <fstream>
#include <stdexcept>

// Implementation for building the state from a .gro file 
void State::parseGRO(const string& groFile) {
    proteinMap.clear();
    residueSet.clear();

    ifstream file(groFile);
    if (!file.is_open()) 
        throw runtime_error("Error opening .gro file: " + groFile);

    // Skip the first two lines
    string line;
    getline(file, line);
    getline(file, line);

    char chain = 'A';
    int chainCount = 0;
    bool isProtein = true; 
    bool hasVel = true;
    shared_ptr<Protein> protein = make_shared<Protein>(chain);

    int resCount = -1; 
    string currentResID;
    shared_ptr<Residue> residue;
    string atomName;

    while (getline(file, line)) {

        // Ignore last line with box dimensions
        if (line.size() < 44) break;

        string resName = trim(line.substr(5, 3)); // NOTE: This trims resName to only 3 characters!
        int tempResNum = stoi(line.substr(0, 5)); // resNum from .pdb not always equal to resCount
        string tempResID = to_string(tempResNum) + resName;

        // Start a new residue if needed
        if (!currentResID.empty() && tempResID != currentResID) {
            residueSet[residue->name].insert(residue);

            if (isProtein) {
                protein->addResidue(residue);

                const bool isAA = Constants::aminoAcids.contains(resName);
                if (!isAA) { // Non-protein residue
                    proteinMap[chain] = protein;
                    isProtein = false;
                    chain = 'Z'; // Z-chain for non-protein residues
                }
                else {  // Still in protein
                    const bool isFinalAtom = Constants::finalAtoms.contains(atomName);
                    if (isFinalAtom) {
                        proteinMap[chain] = protein;
                        if (++chainCount > 25)
                            throw runtime_error("Chain limit exceeded 25 during parseGRO().");
                        chain = Constants::alphabet.at(chainCount);
                        protein = make_shared<Protein>(chain);
                    }
                }
            }

            currentResID.clear();
        }

        // Create residue if none exists for this ID
        if (currentResID.empty()) {
            ++resCount;
            residue = make_shared<Residue>(
                to_string(resCount) + resName, // Correct residue ID
                resName,                       // Residue name
                resCount,                      // Correct residue number
                chain                          // Chain
            );
            currentResID = tempResID;
        }

        // Atom components
        atomName = trim(line.substr(9, 6));

        Vec3D velocity{ 0,0,0 };
        if (hasVel) {
            try {
                velocity = {
                    stod(line.substr(44,8)),
                    stod(line.substr(52,8)),
                    stod(line.substr(60,8))
                };
            }
            catch (...) {
                hasVel = false;
            }
        }

        // Build atom
        residue->addAtom(
            make_shared<Atom>(
                atomName,    
                atomName[0],
                Vec3D(       
                    stod(line.substr(20, 8)),
                    stod(line.substr(28, 8)),
                    stod(line.substr(36, 8))),
                velocity,                      
                residue                        
            )
        );
    }
    file.close();

    // Flush last residue/protein
    if (residue) {
        residueSet[residue->name].insert(residue);
        if (isProtein) protein->addResidue(residue);
    }
    if (isProtein) {
        proteinMap[chain] = protein;
    }

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
        nterm->termini_gpb = Constants::GPBs.at("NTERM");
        auto& cterm = protein->residues.back();
        cterm->termini = true;
        cterm->termini_gpb = Constants::GPBs.at("CTERM");
    }
}