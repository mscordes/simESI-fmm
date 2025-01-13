#include "Core.h"
#include <unordered_map>
#include <iostream>
#include <cmath>      
#include <random>   
#include <vector>     
#include <exception>  
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <filesystem>
#include <memory>

namespace Core {

    // pdb2gmx takes residues types in specific order, exclude termini for now
    const std::vector<std::string> pdb2gmxOrder = { "LYS", "ARG", "ASP", "GLU", "HIS" };

    // Titrable hydrogens of each titrable amino acid (CHARMM36)
    const std::unordered_map<std::string, std::vector<std::string>> hydrogensMap = {
        { "LYS", {"HZ1", "HZ2", "HZ3"} }, { "ARG", {"HH11", "HH12", "HH21", "HH22"} },
        { "ASP", {"HD2"} }, { "GLU", {"HE2"} }, { "HIS", {"HD1", "HE2"} },
        {"NTERM", {"H1", "H2", "H3"}},  {"CTERM", {"HT2"}}
    };

	// Map of pdb2gmx inputs per residue type accounting for termini
	static std::unordered_map<std::string, std::vector<std::string>> createPdb2gmxMap(const Protein& monomer) {
		std::unordered_map<std::string, std::vector<std::string>> pdb2gmxMap;

        // Non-termini 
        pdb2gmxMap.emplace("LYS", std::vector<std::string>{"1", "0"});
        pdb2gmxMap.emplace("ARG", std::vector<std::string>{"1", "0"});
        pdb2gmxMap.emplace("ASP", std::vector<std::string>{"1", "0"});
        pdb2gmxMap.emplace("GLU", std::vector<std::string>{"1", "0"});
        pdb2gmxMap.emplace("HIS", std::vector<std::string>{"2", "0"});

        // Termini, start with N-termini
        const std::string& nterm_res_name = monomer.residues[0].lock()->res_name;
        if (nterm_res_name == "MET") {
            pdb2gmxMap.emplace("NTERM", std::vector<std::string>{"1", "2"});
        }
        else if (nterm_res_name == "PRO") {
            pdb2gmxMap.emplace("NTERM", std::vector<std::string>{"1", "0"});
        }
        else {
            pdb2gmxMap.emplace("NTERM", std::vector<std::string>{"0", "1"});
        }

        // C-termini always same in CHARMM36
        pdb2gmxMap.emplace("CTERM", std::vector<std::string>{"1", "0"});

		return pdb2gmxMap;
	}

	// Protonation probability via Henderson-Hasselbalch equation
    static float protProbability(float pka_val, float pH) {
        return static_cast <float>(1.0 / (1 + std::pow(10, pH - pka_val)));
    }

    static int randomChoice(float prob) {
        std::random_device rd;  
        std::mt19937 gen(rd()); 
        std::discrete_distribution<> d({ prob, 1 - prob });
        return d(gen);  // 0 for protonated, 1 for deprotonated
    }

    // Returns a vector of pdb2gmx inputs for a given residue type within a given protein
    static std::vector<std::string> setProtState(
        const std::string& resType,
        const std::vector<std::shared_ptr<Residue>>& residues,
        const std::unordered_map<std::string, std::vector<std::string>>& pdb2gmxMap,
        const std::unordered_map<int, float>& pkaVals,
        const int& upperbound, const int& lowerbound, const Config& config) {

        std::vector<std::string> residueInputs;

        // Get val for each residue of the inputted type
        for (const auto& residue : residues) {
            if (residue->res_num >= lowerbound && residue->res_num <= upperbound) {
                try {
            		auto it = pkaVals.find(residue->res_num);
                    if (it != pkaVals.end()) {
                        float pka = it->second;
						// FIX - what should the initial pH be?
                        float prob = protProbability(pka, 7.0f);
                        int protState = randomChoice(prob);

                        // Add prot_state to residueInputs or perform other necessary operations
                        auto mapIt = pdb2gmxMap.find(resType);
                        if (mapIt != pdb2gmxMap.end()) {
                            residueInputs.push_back(mapIt->second[protState]);
                        }
                        else {
                            std::cerr << "Residue name not found in pdb2gmx map." << std::endl;
                            std::exit(1);
                        }
                    } else {
                        throw std::runtime_error("Residue number not found in pKa values map.");
                    }
                }
                catch (const std::exception& e) {
                    std::cerr << "\n\n\nERROR: PROPKA .pka file was generated, but is likely missing a residue.\n"
                        "A common issue is that the C-termini for each monomer must have both oxygens, "
                        "and the final oxygen must have the atom name 'OXT', not 'OT1' or 'OT2'."
                        << std::endl;
                    std::cerr << "Exception: " << e.what() << std::endl;
					std::cerr << "Residue number: " << residue->res_num << std::endl;
                    std::exit(1);
                }
            }
        }
        return residueInputs;

    }

	/* Creates vectors of ints to feed into pdb2gmx to generate.top files with protonation states
	   corresponding to the pKa values of the residues in the protein at pH 5.26 (+ESI) or 8.75 (-ESI) */
    std::unordered_map<std::string, std::vector<std::string>> setProtStates(const CoordInfo& coordInfo, const std::unordered_map<int, float>& pkaVals,
        const Config& config) {

        std::unordered_map<std::string, std::vector<std::string>> pdb2gmxInputs;

        // Each monomer will have its own vector of inputs
        for (const auto& monomer : coordInfo.proteins) {
			std::vector<std::string> monomerInputs;

            // Dict of pdb2gmx prot state inputs corrected for termini
			std::unordered_map<std::string, std::vector<std::string>> pdb2gmxMap = createPdb2gmxMap(*monomer);

            // Bounds to determine if residues within this monomer
            int upperBound = monomer->residues.back().lock()->res_num;
            int lowerBound = monomer->residues[0].lock()->res_num;

            // pdb2gmx takes residues types in specific order, exclude termini for now
			for (const auto& resType : pdb2gmxOrder) {
                if (coordInfo.residueMap.find(resType) != coordInfo.residueMap.end()) {
                    std::vector<std::string> resInputs = setProtState(resType, coordInfo.residueMap.at(resType), 
                        pdb2gmxMap, pkaVals, upperBound, lowerBound, config);
                    monomerInputs.insert(monomerInputs.end(), resInputs.begin(), resInputs.end());
                }
            }

            // Add termini
            std::vector<std::shared_ptr<Residue>> nterm = { monomer->residues[0].lock()};
            std::vector<std::string> ntermInput = setProtState("NTERM", nterm, pdb2gmxMap, pkaVals, upperBound, 
                lowerBound, config);
            monomerInputs.insert(monomerInputs.end(), ntermInput.begin(), ntermInput.end());

            std::vector<std::shared_ptr<Residue>> cterm = { monomer->residues.back().lock()};
            std::vector<std::string> ctermInput = setProtState("CTERM", cterm, pdb2gmxMap, pkaVals, upperBound, 
                lowerBound, config);
            monomerInputs.insert(monomerInputs.end(), ctermInput.begin(), ctermInput.end());

            // Add monomer
			pdb2gmxInputs[monomer->chain] = monomerInputs;
        }

		return pdb2gmxInputs;
    }

	// Determine protonation state of a resiude type, and return the corresponding pdb2gmx inputs
    static std::vector<std::string> getProtState(const std::string& resType, const std::vector<std::shared_ptr<Residue>>& residues,
        const std::unordered_map<std::string, std::vector<std::string>>& pdb2gmxMap, const int& upperbound, const int& lowerbound,
        const std::vector<std::string>& hydrogens) {

        std::vector<std::string> residueInputs;

        // Get val for each residue of the inputted type
        for (const auto& residue : residues) {
            if (residue->res_num >= lowerbound && residue->res_num <= upperbound) {
				std::vector<std::string> hydrogensPresent;
                
				for (const auto& weak_atom : residue->atoms) {
                    auto atom = weak_atom.lock();
					if (std::find(hydrogens.begin(), hydrogens.end(), atom->atom_name) != hydrogens.end()) {
						hydrogensPresent.push_back(atom->atom_name);
					}
				}

				if (hydrogensPresent.size() == hydrogens.size()) {
					residueInputs.push_back(pdb2gmxMap.at(resType)[0]);
				}
				else {                    
                    // HIS can have 3 protonation states, so parse here
                    if (residue->res_name == "HIS") {
                        for (const auto& weak_atom : residue->atoms) {
                            auto atom = weak_atom.lock();
                            if (atom->atom_name == "HD1") {
                                residueInputs.push_back("0"); // HISD
                                break;
                            }
                            else if (atom->atom_name == "HE2") {
                                residueInputs.push_back("1"); // HISE
                                break;
                            }
                        }
                    }
					else {
						residueInputs.push_back(pdb2gmxMap.at(resType)[1]);
					}
				}
            }
        }
        
        return residueInputs;
    }

    // Create map of pdb2gmx inputs of the protonation states of the inputted protein atoms 
    std::unordered_map<std::string, std::vector<std::string>> getProtStates(const CoordInfo& coordInfo) {
        std::unordered_map<std::string, std::vector<std::string>> pdb2gmxInputs;

        // Each monomer will have its own vector of inputs
        for (const auto& monomer : coordInfo.proteins) {
            std::vector<std::string> monomerInputs;

            // Dict of pdb2gmx prot state inputs corrected for termini
            std::unordered_map<std::string, std::vector<std::string>> pdb2gmxMap = createPdb2gmxMap(*monomer);

            // Bounds to determine if residues within this monomer
            int upperBound = monomer->residues.back().lock()->res_num;
            int lowerBound = monomer->residues[0].lock()->res_num;

            // pdb2gmx takes residues types in specific order, exclude termini for now
            for (const auto& resType : pdb2gmxOrder) {
                if (coordInfo.residueMap.find(resType) != coordInfo.residueMap.end()) {
                    std::vector<std::string> resInputs = getProtState(resType, coordInfo.residueMap.at(resType), pdb2gmxMap, 
                        upperBound, lowerBound, hydrogensMap.at(resType));
                    monomerInputs.insert(monomerInputs.end(), resInputs.begin(), resInputs.end());
                }
            }

            // Add termini
            std::vector<std::shared_ptr<Residue>> nterm = { monomer->residues[0].lock() };
            std::vector<std::string> ntermInput = getProtState("NTERM", nterm, pdb2gmxMap, upperBound, lowerBound, hydrogensMap.at("NTERM"));
            monomerInputs.insert(monomerInputs.end(), ntermInput.begin(), ntermInput.end());

            std::vector<std::shared_ptr<Residue>> cterm = { monomer->residues.back().lock() };
            std::vector<std::string> ctermInput = getProtState("CTERM", cterm, pdb2gmxMap, upperBound, lowerBound, hydrogensMap.at("CTERM"));
            monomerInputs.insert(monomerInputs.end(), ctermInput.begin(), ctermInput.end());

            // Add monomer
            pdb2gmxInputs[monomer->chain] = monomerInputs;
        }

        return pdb2gmxInputs;
    }

    // Make .top file with proper .itp info and residue ordering/numbering
    void createTOP(const std::string& fname, const CoordInfo& coordInfo, const std::vector<std::string>& topOrder) {
        std::ofstream file(fname);

        if (!file.is_open()) {
            std::exit(1);
        }

        // Begin writing .top, start with .itp files
        file << "#include \"charmm36.ff/forcefield.itp\"\n"; // FF info
        // Protein .itp files
		for (const auto& monomer : coordInfo.proteins) {
			file << "#include \"" << monomer->chain << ".itp\"\n";
		}
		file << "#include \"o2.itp\"\n";            // O2 gas
        file << "#include \"n2.itp\"\n";            // N2 gas
		file << "#include \"nh4.itp\"\n";           // Ammonium
        file << "#include \"nh3.itp\"\n";           // Ammonia
        file << "#include \"aceh.itp\"\n";          // Acetic acid
        file << "#include \"ace.itp\"\n";           // Acetate
        file << "#include \"hydroxide.itp\"\n";     // Hydroxide ions
		file << "#include \"hydronium.itp\"\n";     // Hydronium ions
		file << "#include \"tip4p_2005.itp\"\n\n";  // TIP4P water model

        // System and molecule information
        file << "[ system ]\n";
        file << "simESI System\n\n";
		file << "[ molecules ]\n";  

        // Protein(s)
        for (const auto& monomer : coordInfo.proteins) {
            file << monomer->chain << "                   " << "1\n";
        }

		// Other molecules
		for (const auto& resname : topOrder) {
			if (coordInfo.numResidues.find(resname) != coordInfo.numResidues.end()) {
				file << resname << "                 " << coordInfo.numResidues.at(resname) << "\n";
			}
		}
         
		file.close();
    }

	// Convert .top file of a given protein chain to an .itp file
    static void convert_top2itp(std::string chain) {
        std::vector<int> indices_toDelete;
        std::vector<std::string> top;

        // Open the .top file for reading
		std::string topFile = chain + ".top";
        std::ifstream infile(topFile);
        std::string line;

        if (!infile.is_open()) {
            std::cerr << "Error processing .top file for protein chain" << chain << 
                ".\nLikely caused by failure of pdb2gmx due to improper input structure." << std::endl;
            std::exit(1);
        }

        // Read all lines from the .top file
        while (std::getline(infile, line)) {
            top.push_back(line);
        }
        infile.close();

        // Find the lines to delete from the bottom of the file
        for (size_t indx = top.size(); indx-- > 0;) {
            if (top[indx].find("; Include Position restraint file") == std::string::npos) {
                indices_toDelete.push_back(static_cast<int>(indx));
            }
            else {
                indices_toDelete.push_back(static_cast<int>(indx));
                break;
            }
        }

        // Find the lines to delete from the top of the file
        for (size_t indx = 0; indx < top.size(); ++indx) {
            if (top[indx].find("[ moleculetype ]") == std::string::npos) {
                indices_toDelete.push_back(static_cast<int>(indx));
            }
            else {
                break;
            }
        }

        // Sort the indices in reverse order to delete lines
        std::sort(indices_toDelete.rbegin(), indices_toDelete.rend());

        // Delete the lines from the file
        for (int indx : indices_toDelete) {
            top.erase(top.begin() + indx);
        }

        // Modify the line with the protein information
        top[2] = chain + "          3";

        // Write the modified content to a new .itp file
        std::ofstream outfile(chain + ".itp");
        for (const auto& line : top) {
            outfile << line << "\n";
        }
        outfile.close();

        // Delete the obsolete .top file
        if (!std::filesystem::remove(topFile)) {
            std::cout << "File " << topFile << " not found.\n";
            std::exit(1);
        }
    }

    // Creates .gro, .top and .itp file(s) given protein protonation states as defined by input
    void writeTOP(const std::string& topName, const CoordInfo& coordInfo, const std::unordered_map<std::string,
        std::vector<std::string>>&pdb2gmxMap, const bool keepGro, const Config& config, const std::vector<std::string>& topOrder) {

        // Create unique .itp files for each monomer
        for (const auto& monomer : coordInfo.proteins) {

            // Make monomer specific .gro file
            const std::string new_groName = "top_" + monomer->chain + ".gro";
			const std::string pre_groName = "pre_" + monomer->chain + ".gro";
			writeGRO(pre_groName, coordInfo.proteinAtoms.at(monomer->chain), coordInfo.box_vectors);

			// Make monomer specific .top file
            {
                const std::string new_topName = monomer->chain + ".top";
                auto it = pdb2gmxMap.find(monomer->chain);
                if (it == pdb2gmxMap.end()) {
                    std::cerr << "Chain not found in pdb2gmxMap: " << monomer->chain << std::endl;
                    std::exit(1);
                }
                const std::vector<std::string>& inputs = it->second;

                std::ostringstream command;
                deleteFile(std::filesystem::path(new_topName));
                command << config.gmx_env << " pdb2gmx -f " << pre_groName << " -o " << new_groName << " -p " << new_topName
                    << " -ff charmm36 -water none -lys -arg -glu -asp -his -ter -ignh";
                auto_gmx_input(command.str(), inputs);

				// Check that .top file was created
				if (!std::filesystem::exists(new_topName)) {
					std::cerr << "Error: .top file not created for protein chain " << monomer->chain << std::endl;
					std::exit(1);
				}

				// Delete the unneccesary posre and .gro files (if not keeping)
                if (!keepGro) {
                    deleteFile(std::filesystem::path(new_groName));
                }
                deleteFile(std::filesystem::path(pre_groName));
                deleteFile(std::filesystem::path("posre.itp"));
            }

            // Convert .top to .itp
			convert_top2itp(monomer->chain);
        }

        // If keeping .gro, need to stitch monomers together
        if (keepGro) {
            std::vector<std::shared_ptr<Atom>> atoms;
            for (const auto& monomer : coordInfo.proteins) {
                std::string monomerGro = "top_" + monomer->chain + ".gro";
                std::vector<std::shared_ptr<Atom>> monomerAtoms = readGROatoms(monomerGro); //.gro files
                atoms.insert(atoms.end(), monomerAtoms.begin(), monomerAtoms.end());

                // Delete the now obsolete monomer .gro file
                if (!std::filesystem::remove(monomerGro)) {
                    std::cout << "File " << monomerGro << " not found.\n";
                    std::exit(1);
                }
            }

            // Final complete .gro
            writeGRO("system.gro", atoms, coordInfo.box_vectors);
        }

        // Create the master .top file
		createTOP(topName, coordInfo, topOrder);
	}

    // Reorder atoms to correspond to the order in the .top file, rebuilds coordInfo in place
    void reorderAtoms(CoordInfo& coordInfo, const std::vector<std::string>& topOrder) {

        // Bin (updated) atoms into Residues
        auto residues = binAtoms(coordInfo.atoms);

        // Get Residue Map
        auto residueMap = getResidueMap(residues);

        // Bin Residues into proteins (more accurately chains)
        auto proteins = binResidues(residues);

        // Get protein atoms, specifically
        auto proteinAtoms = getProteinAtoms(proteins);

        std::vector<std::shared_ptr<Atom>> newAtoms;
        int resCount = -1;
		int atomCount = 0;

        // Initialize res_id
        std::string temp_id = "foo";

        // Proteins go first
        for (const auto& pair : proteinAtoms) {
			for (const auto& atom : pair.second) {

				// Determine if new residue
				if (atom->res_id != temp_id) {
					resCount++;
					temp_id = atom->res_id;
				}

				std::shared_ptr<Atom> newAtom = std::make_shared<Atom>(*atom);
				newAtom->atom_num = atomCount;
				newAtom->res_num = resCount;
				newAtom->res_id = std::to_string(resCount) + newAtom->res_name;
                newAtoms.push_back(newAtom);
				atomCount++;
            }
        }

        // Non-protein residues
        for (const auto& resType : topOrder) {
            auto it = residueMap.find(resType);

            if (it != residueMap.end() && it->second.size() > 0) {
            	for (const auto& residue : it->second) {
                    for (auto& weak_atom : residue->atoms) {
                        auto atom = weak_atom.lock();

                        // Determine if new residue
                        if (atom->res_id != temp_id) {
                            resCount++;
                            temp_id = atom->res_id;
                        }

                        atom->atom_num = atomCount;
                        atom->res_num = resCount;
                        atom->res_id = std::to_string(resCount) + atom->res_name;
                        newAtoms.push_back(atom);
                        atomCount++;
                    }
            	}
            }
        }

        // Rebuild coordInfo with now properly ordered atoms
        coordInfo = parseAtoms(newAtoms, coordInfo.box_vectors[0]);
    }

	// Takes input .gro file and outputs .top file preserving protonation states
    void topFromCoord(const std::string& groFname, const std::string& topFname, const float& boxSize, 
        const Config& config, const std::vector<std::string>& topOrder) {
		
        // Get protonation states from .gro file
		CoordInfo coordInfo = buildCoordInfo(groFname, boxSize);
		const std::unordered_map<std::string, std::vector<std::string>> pdb2gmxMap = getProtStates(coordInfo);

		// Write the .top file
		writeTOP(topFname, coordInfo, pdb2gmxMap, false, config, topOrder);
    }
}