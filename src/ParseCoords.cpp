#include "Core.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <memory>

namespace Core {

    static std::string removeSpaces(std::string input) {
        input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
        return input;
    }

    // Returns a list vector of Atom class object pointers parsed from an input .pdb
    static std::vector<std::shared_ptr<Atom>> readPDBatoms(const std::string& pdb) {
        std::vector<std::shared_ptr<Atom>> atoms; 
        std::ifstream file(pdb);

        if (!file.is_open()) {
            std::cerr << "Error opening .pdb file: " << pdb << std::endl;
            std::exit(1);
        }

        int atomCount = 0;
		int resCount = -1;
		int chainCount = 0;
		std::string temp_res_id = "0";
		std::string temp_chain = "foo";
        std::string line;
        while (std::getline(file, line)) {
            // Check if the line starts with "ATOM"
            if (line.substr(0, 4) == "ATOM") {
                // Determine if new residue for residue numbering
                std::shared_ptr<Residue> parent = nullptr;                       // Parent Residue object, leave for now
                std::string res_name = removeSpaces(line.substr(17, 3));         // Residue Name
                int res_num = std::stoi(line.substr(22, 4));                     // Residue number from .pdb
                std::string test_res_id = std::to_string(res_num) + res_name;    // Residue ID from .pdb
				if (test_res_id != temp_res_id) {
					resCount += 1;
					temp_res_id = test_res_id;
				}
                res_num = resCount;                                         // Corrected residue number
                std::string res_id = std::to_string(res_num) + res_name;    // Corrected residue ID  

                std::string atom_name = removeSpaces(line.substr(12, 4));   // Atom name 
                int atom_num = atomCount;                                   // Atom number

                // Chain, convert from alphabetical to numerical
				std::string chain = std::to_string(chainCount); 
				std::string test_chain = removeSpaces(line.substr(21, 1)); // Chain ID from .pdb
                if (test_chain != temp_chain) {
					chainCount++;
					temp_chain = test_chain;
                    chain = std::to_string(chainCount);
                }

                // Extract coordinates (in nm)
                std::array<float, 3> coord = {
                 std::stof(line.substr(30, 8)) / 10,
                 std::stof(line.substr(38, 8)) / 10,
                 std::stof(line.substr(46, 8)) / 10
                };

                // Initialize velocity as zero 
                std::array<float, 3> velocity = { 0.0f, 0.0f, 0.0f }; 

                std::string element = line.substr(76, 2); // Element

                // Create an Atom object and add it to the list
				atoms.push_back(std::make_shared<Atom>(parent, res_id, res_num, res_name, atom_name, atom_num, coord, velocity, element, chain));

                // Count up iterator
                ++atomCount;
            }
        }

        file.close(); 
        return atoms;  
    }

    // Returns a list vector of Atom class object pointers parsed from an input .gro
    std::vector<std::shared_ptr<Atom>> readGROatoms(const std::string& gro) {
        std::vector<std::shared_ptr<Atom>> atoms;
        std::ifstream file(gro);

        if (!file.is_open()) {
            std::cerr << "Error opening .gro file: " << gro << std::endl;
            std::exit(1);
        }

        std::string line;
        std::vector<std::string> lines;

        // Skip the first two lines
        std::getline(file, line);
        std::getline(file, line); 
        
		// Use num atoms stored in .gro to reserve space for atoms
		int numAtoms;
        try {
            numAtoms = std::stoi(line);
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "Improperly formatted .gro file." << e.what() << std::endl;
            std::exit(1);
        }
        catch (const std::out_of_range& e) {
            std::cerr << "Improperly formatted .gro file." << e.what() << std::endl;
            std::exit(1);
        }
		atoms.reserve(numAtoms);

        // Read the rest of the lines
        while (std::getline(file, line)) {
            lines.push_back(line);
        }

        // Determine if .gro has velocities
		bool hasVel = splitLine(lines[0]).size() == 9;

        // Process .gro, but skip last line
        int atomCount = 0;
		int resCount = -1;
        int chainCount = 1;
        std::string chain(std::to_string(chainCount)); // Chain
        std::string temp_res_id = "0";
		bool end_prot = false;
        for (size_t i = 0; i < lines.size() - 1; ++i) {
            // Determine if new residue for residue numbering
            std::string res_name = removeSpaces(lines[i].substr(5, 3));     // Residue Name (only take 3 letter code) 
            int res_num = std::stoi(lines[i].substr(0, 5));                 // Residue number from .gro
            std::string test_res_id = std::to_string(res_num) + res_name;   // Residue ID from .gro
            if (test_res_id != temp_res_id) {
                resCount += 1;
                temp_res_id = test_res_id;
            }
			res_num = resCount;                                             // Corrected residue number
            std::string res_id = std::to_string(res_num) + res_name;        // Corrected residue ID  

            std::string atom_name = removeSpaces(lines[i].substr(9, 6));    // Atom name 
            int atom_num = atomCount;                                       // Atom number
            std::string element(1, atom_name[0]);                           // Element

            // Update chain to account for complexes/non-protein residues
            if (end_prot) {

                // Non-protein atoms
                if (res_name == "SOL" || res_name == "NNN" || res_name == "HHO" ||
                    res_name == "OHX" || res_name == "ATX" || res_name == "AHX" ||
                    res_name == "NXX" || res_name == "NXH" || res_name == "OOO") {
					chain = "Z"; // Z-chain for non-protein
                }
                end_prot = false;
            }

            // Find end of protein since .gro files don't hold chain info
            if (atom_name == "OXT" || atom_name == "OT2" || atom_name == "HT2") {
                // Determine if actual end of protein
                if (lines[i + 1].find("HT2") == std::string::npos) {
                    end_prot = true;
                    chainCount++;
					chain = std::to_string(chainCount);
                }
            }

            // Extract coordinates
            std::array<float, 3> coord = {
                std::stof(lines[i].substr(20, 8)),
                std::stof(lines[i].substr(28, 8)),
                std::stof(lines[i].substr(36, 8))
            };
            
            // Extract velocity 
            std::array<float, 3> velocity = { 0.0f, 0.0f, 0.0f };
            if (hasVel) {
                velocity = {
                    std::stof(lines[i].substr(44, 8)),
                    std::stof(lines[i].substr(52, 8)),
                    std::stof(lines[i].substr(60, 8))
                };

                // Velocity check
				if (velocity[0] > 10.0f || velocity[1] > 10.0f || velocity[2] > 10.0f) {
                    if (element != "M") {
                        velocity = sampleMaxwell(1, 300.0f, element)[0];
                    }
				}
            }

            // Create an Atom object and add it to the list
			atoms.push_back(std::make_shared<Atom>(nullptr, res_id, res_num, res_name, atom_name, atom_num, coord, velocity, element, chain));
            
            // Count up iterator
            ++atomCount;
        }

        file.close();
        return atoms;
    }

    // Bins Atom objects into Residues
    std::vector<std::shared_ptr<Residue>> binAtoms(std::vector<std::shared_ptr<Atom>>& atoms) {
        std::vector<std::shared_ptr<Residue>> residues;
        if (!atoms.empty())
		    residues.reserve(atoms.back()->res_num + 1);

        // Initialize the first residue
        std::shared_ptr<Residue> residue = std::make_shared<Residue>(atoms[0]->res_id, 0, atoms[0]->res_name, atoms[0]->chain);
        std::string temp_id = atoms[0]->res_id;

        int count = 1;
        for (auto& atom : atoms) {
            // If res_id changes, push current residue and start a new one
            if (atom->res_id != temp_id) {
                residues.push_back(residue);

                // Update each atom in residues partent to point to the residue
				std::weak_ptr<Residue> weak_residue = residue;
				for (const auto& res_atom : residue->atoms) {
                    res_atom.lock()->parent = weak_residue;
				}

                // Start new residue
                residue = std::make_shared<Residue>(atom->res_id, count, atom->res_name, atom->chain);
                temp_id = atom->res_id;
                count += 1;
            }
            residue->addAtom(atom);
        }

        // Push the last residue into the vector and update child atom objects
        residues.push_back(residue);
        std::weak_ptr<Residue> weak_residue = residue;
        for (const auto& res_atom : residue->atoms) {
            res_atom.lock()->parent = weak_residue;
        }

        return residues;
    }

    // Bins Residue objects into Proteins
    std::vector<std::shared_ptr<Protein>> binResidues(const std::vector<std::shared_ptr<Residue>>& residues) {
        std::vector<std::shared_ptr<Protein>> proteins;

        // Initialize the first protein
        std::shared_ptr<Protein> protein = std::make_shared<Protein>(residues[0]->chain);
        std::string temp_chain = residues[0]->chain;

        for (const auto& residue : residues) {
            // Z-chain is non-protein so ignore 
			if (residue->chain == "Z") {
				break;
			}

            // If chain changes, push current protein and start a new one
            else if (residue->chain != temp_chain) {
                proteins.push_back(protein);
                protein = std::make_shared<Protein>(residue->chain);
                temp_chain = residue->chain;
            }
            protein->addResidue(residue);
        }

        // Push the last protein into the vector
        proteins.push_back(protein);

        return proteins;
    }

    // Returns a map of protein monomers where val corresponds to atoms in that monomer
    std::map<std::string, std::vector<std::shared_ptr<Atom>>> getProteinAtoms(const std::vector<std::shared_ptr<Protein>>& proteins) {
        std::map<std::string, std::vector<std::shared_ptr<Atom>>> proteinAtoms;

        // Parse monomers
        for (const auto& monomer : proteins) {
            std::vector<std::shared_ptr<Atom>> monomerAtoms;

            // Concatenate atoms from each residue
            for (const auto& residue : monomer->residues) {
                for (const auto& atom : residue.lock()->atoms) {
                    monomerAtoms.push_back(atom.lock());
                }
            }
            
            // Add monomer atoms to map
            proteinAtoms[monomer->chain] = std::move(monomerAtoms); 
        }

        return proteinAtoms;
    }

    // Returns array of atomic coordinated given a vector of Atom objects
    std::vector<std::array<float, 3>> extractCoordinates(const std::vector<std::shared_ptr<Atom>>& atoms) {
        std::vector<std::array<float, 3>> coordinates;
        coordinates.reserve(atoms.size()); 

        for (const auto& atom : atoms) {
            coordinates.push_back(atom->coord); 
        }
        
        return coordinates;
    }

    // Returns a map of protein monomers where val corresponds to coords of atoms in that monomer
    std::vector<std::array<float, 3>> getProteinCoords(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms) {
       std::vector<std::array<float, 3>> proteinCoords;

        for (const auto& pair : proteinAtoms) {
			for (const auto& atom : pair.second) {
				proteinCoords.push_back(atom->coord);
			}
        }

        return proteinCoords;
    }

    /* Returns a Residue Map, where each key is a residue name (like Arg), and val = list of Resiude objects
    of that type given an input list of all residues in the system */
    std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>> 
        getResidueMap(const std::vector<std::shared_ptr<Residue>>& residues) {

        std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>> residueMap;

        // Vector of the all important residues
        const std::vector<std::string> titResnames = {
            "LYS", "ARG", "GLU", "ASP", "HIS", "SOL", "HHO", "OHX",
            "ATX", "AHX", "NXX", "NXH", "NNN", "OOO"
        };
        residueMap.reserve(titResnames.size() + residues.size());

        // Group Residues based on residue name
        for (const auto& residue : residues) {
            residueMap[residue->res_name].push_back(residue);
        }

        // Add in residue types not present
        for (const auto& res_name : titResnames) {
            if (residueMap.find(res_name) == residueMap.end()) {
                residueMap[res_name] = {};
            }
        }

        return residueMap;
    }

    // Returns a map with the number of each type of residue in the system given a list of residues
    std::unordered_map<std::string, int> getNumResidues(
        const std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>>& residueMap) {

        std::unordered_map<std::string, int> numResidues;
        numResidues.reserve(residueMap.size());

        for (const auto& pair : residueMap) {
            numResidues.emplace(pair.first, static_cast<int>(pair.second.size()));
        }

        return numResidues;
    }

    // Get water O coordinates
	static std::vector<std::array<float, 3>> getWaterOCoords(const std::vector<std::shared_ptr<Residue>>& waters) {
		std::vector<std::array<float, 3>> waterOCoords;
		waterOCoords.reserve(waters.size());

		for (const auto& water : waters) {
            waterOCoords.push_back(water->atoms[0].lock()->coord);
		}

		return waterOCoords;
	}

    // Get water H coordinates
    static std::vector<std::array<float, 3>> getWaterHCoords(const std::vector<std::shared_ptr<Residue>>& waters) {
        std::vector<std::array<float, 3>> waterHCoords;
        waterHCoords.reserve(waters.size());

        for (const auto& water : waters) {
            waterHCoords.push_back(water->atoms[1].lock()->coord); //HW2
            waterHCoords.push_back(water->atoms[2].lock()->coord); //HW3
        }

        return waterHCoords;
    }

    // Extracts info from vector of atoms and outputs CoordInfo object
    CoordInfo parseAtoms(std::vector<std::shared_ptr<Atom>> atoms, const float& boxSize) {
        // Bin Atoms into Residues
        auto residues = binAtoms(atoms);

        // Bin Residues into Proteins 
        auto proteins = binResidues(residues);

        // Extract coordinates of all atoms
        auto coords = extractCoordinates(atoms);

        // Get atoms specific to protein(s)
        auto proteinAtoms = getProteinAtoms(proteins);

        // Get atom coordinates specific to protein(s)
        auto proteinCoords = getProteinCoords(proteinAtoms);

        // Get simulation box vectors
        std::array<float, 3> boxVectors = { boxSize, boxSize, boxSize };

        // Get Residue dict
        auto residueMap = getResidueMap(residues);

        // Get number of each residue type
        auto numResidues = getNumResidues(residueMap);

        // Get water O coordinates
        auto waterOCoords = getWaterOCoords(residueMap["SOL"]);

        // Get water H coordinates
        auto waterHCoords = getWaterHCoords(residueMap["SOL"]);

        // Build CoordInfo object
        return CoordInfo(
            std::move(atoms),
            std::move(residues),
            std::move(proteins),
            std::move(proteinAtoms),
            std::move(coords),
            std::move(proteinCoords),
            std::move(waterOCoords),
            std::move(waterHCoords),
            boxVectors,
            std::move(residueMap),
            std::move(numResidues)
        );
    }


    /* Combines previous functions and outputs a useful CoordInfo object with combined information
    for a particular coordinate file */
    CoordInfo buildCoordInfo(const std::string& file, const float& boxSize) {
        // Determine file type and unpack into Atoms
        std::string ftype; 
		if (file.length() < 5) {
			std::cerr << "Invalid inputted .pdb." << std::endl;
			std::exit(1);
		}
		else {
			ftype = file.substr(file.length() - 4);
		}
        std::vector<std::shared_ptr<Atom>> atoms;
        if (ftype == ".gro") {
            atoms = readGROatoms(file); //.gro files
        }
        else if (ftype == ".pdb") {
            atoms = readPDBatoms(file); //.pdb files
        }
		else {
			std::cerr << "Invalid inputted file type." << std::endl;
			std::exit(1);
		}

		// Parse Atoms into CoordInfo object
		CoordInfo coordInfo = parseAtoms(atoms, boxSize);

		return coordInfo;
    }
}