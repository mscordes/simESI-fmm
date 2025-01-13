#include "Core.h"
#include <iostream>
#include <sstream>
#include <cstring> 
#include <unordered_map>

namespace Core {

    // Can still pull other non-pKa information so delete fluff here
    static bool isInteger(const std::string& s) {
        std::stringstream ss(s);
        int x;
        return (ss >> x) && (ss.eof());  // Returns true if the entire string is converted to an integer
    }

    // Generate a .pka file and return a map of pKa values corresponding to chargeable amino acids
    std::unordered_map<int, float> getPkaVals(const CoordInfo& coordInfo, const std::string& pdb) {
        // Generate .pka via PROPKA3
        std::string command = "propka3 " + pdb + ".pdb";
        std::vector<std::string>inputs = {};
        auto_gmx_input(command, inputs);

        // Open .pka file
        std::string pkaFile = std::string(pdb) + ".pka";
		if (!std::filesystem::exists(pkaFile)) {
			std::cerr << "pKa file not found. Likely do not have PROPKA3 installed and callable from command line." << std::endl;
			exit(1);
		}
        std::vector<std::string> fileLines = Core::readFile(pkaFile);
        std::vector<std::vector<std::string>> splitFileLines = splitLines(fileLines);

        // Adjust residue numbering if complex
        int protCount = 0;
		std::unordered_map<std::string, int> resnumMap;
        for (auto& pair : coordInfo.proteinAtoms) {
            
			// Only 26 letters in alphabet so max 26 chains, this should be fixed later
            int chainNum = std::stoi(pair.first);
			std::cout << "Chain number: " << chainNum << std::endl; 
			if (chainNum > 26) {
				std::cerr << "ERROR: More than 26 chains detected. PROPKA3 (maybe?) only supports 26 chains." << std::endl;
				exit(1);
			}
            std::string chain(1, static_cast<char>(chainNum - 1) + 'A');
			resnumMap[chain] = pair.second[0]->res_num;
            protCount++;
        }

        // Create map
        std::unordered_map<int, float> pkaMap;
        for (const auto& lineWords : splitFileLines) {
            if (lineWords.size() == 5 && isInteger(lineWords[1])) {
				std::string chain = lineWords[2];
                int resnum = std::stoi(lineWords[1]) - 1 + resnumMap[chain]; // Account for 0-based indexing 
                float pka = std::stof(lineWords[3]);  
                pkaMap[resnum] = pka;  
            }
        }

        return pkaMap;
    }

	// Updates pKa values based on current conformation
    void updatePkavals(std::unordered_map<int, float>& pkaMap, CoordInfo& coordInfo) {

        // Final oxygens must be of the type OXT or PROPKA throws an error
        std::unordered_map<std::string, std::string> tempNames = {};
        for (const auto& monomer : coordInfo.proteinAtoms) {
            std::shared_ptr<Atom> lastAtom = monomer.second.back();
            std::shared_ptr<Atom> secondLastAtom = monomer.second.end()[-2];

            if (lastAtom->element == "O") {
                tempNames[monomer.first] = lastAtom->atom_name;
                lastAtom->atom_name = "OXT";
            }
			else if (secondLastAtom->element == "O") {
                tempNames[monomer.first] = secondLastAtom->atom_name;
				secondLastAtom->atom_name = "OXT";
			}
            else {
				std::cerr << "ERROR: Could not find final oxygen in protein for PROPKA to rename." << std::endl;
				exit(1);
            }
        }

        // Write new .pdb file with corrected naming and update pKa values
        writePDB("pka.pdb", coordInfo.proteinAtoms, coordInfo.box_vectors);
        pkaMap = getPkaVals(coordInfo, "pka");

        // Reset final atom names
        for (const auto& monomer : coordInfo.proteinAtoms) {
            std::shared_ptr<Atom> lastAtom = monomer.second.back();
            std::shared_ptr<Atom> secondLastAtom = monomer.second.end()[-2];

            if (lastAtom->element == "O") {
                lastAtom->atom_name = tempNames[monomer.first];
            }
            else if (secondLastAtom->element == "O") {
                secondLastAtom->atom_name = tempNames[monomer.first];
            }
            else {
                std::cerr << "ERROR: Could not find final oxygen in protein to revert PROPKA naming." << std::endl;
                exit(1);
            }
        }
    }

}