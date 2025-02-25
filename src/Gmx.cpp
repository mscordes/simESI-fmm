#include "Core.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <algorithm>

#ifdef _WIN32
#include <direct.h>
#define popen _popen
#define pclose _pclose
#endif

namespace Core {

    // Automate command line inputs to command line calls, can also call w/out inputs
    void auto_gmx_input(const std::string& command, const std::vector<std::string>& inputs) {
        FILE* pipe = popen(command.c_str(), "w");
        if (!pipe) {
            std::cerr << "Failed to run command" << std::endl;
            return;
        }

        // If there are inputs, write them to the process's stdin
        for (const auto& input : inputs) {
            fputs((input + "\n").c_str(), pipe);
        }

        // Close the pipe
        if (pclose(pipe) == -1) {
            std::cerr << "Failed to close command pipe: " << command << std::endl;
        }
    }

    // Call mdrun, modify for GPU acceleration or HPC mode dependent on user args
	void runMD(const std::string& fname, const Config& config) {
        std::ostringstream command;
        if (config.hpc == "yes") {
			command << config.gmx_env << " mdrun -ntmpi 1 -nb gpu -v -deffnm " + fname;
		}
		else if (config.gpu == "yes") {
			command << config.gmx_env << " mdrun -nb gpu -v -deffnm " + fname;
		}
        else {
            command << config.gmx_env << " mdrun -v -deffnm " + fname;
        }

        std::vector<std::string> inputs = {};
        auto_gmx_input(command.str(), inputs);
	}

    // Write an index group in an .ndx file, returns bool whether grp is actually present or not
	static bool writeNDXGroup(std::ofstream& file, const std::vector<std::string>& topOrder, const std::string& grp_name, 
        const std::unordered_map<std::string, std::vector<int>>& indices, 
        const std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>>& residueMap, 
        bool addProtein, const std::vector<int>& proteinAtomIndices) {

        // First determine if group present at all
        bool writeGrp = false;
        if (addProtein) {
            writeGrp = true;
        }
        for (const auto& resType : topOrder) {
            if (!indices.at(resType).empty()) {
                writeGrp = true;
                break;
            }
        }
        if (writeGrp) {
            file << "[ " << grp_name << " ]\n";
            int count = 1;

			// If adding protein atoms (for protein + droplet group)
            if (addProtein) {
                int count = 1;
                for (const auto& index : proteinAtomIndices) {
                    file << std::right << std::setw(4) << index << " ";
                    count++;
                    if (count % 16 == 0) {
                        file << '\n';
                        count = 1;
                    }
                }
            }

            // Non-protein atoms
            for (const auto& resType : topOrder) {
                if (indices.at(resType).empty()) {
                    continue;
                }
                for (const auto& index : indices.at(resType)) {
                    for (const auto& atom : residueMap.at(resType)[index]->atoms) {
						auto atom_num = atom.lock()->atom_num + 1;
                        file << std::right << std::setw(4) << atom_num << " ";
                        count++;
                        if (count % 16 == 0) {
                            file << '\n';
                            count = 1;
                        }
                    }
                }
            }
            if (count != 1) {
                file << '\n';
            }
        }

		return writeGrp;
	}   

	// Create .ndx and modify .mdp file for COM translational/rotational removal, T-coupling, FMM, and non-bonded cutoffs
    void modifyMDPgrps(const std::string& ndx_fname, const std::string& mdp_fname, const Config& config, 
        const std::vector<std::string>& topOrder, CoordInfo& coordInfo, const float& temperature) {

		// Get protein atom indices first, also store protein C coords for distance calculations later
		std::vector<int> proteinAtomIndices;
		std::vector<std::array<float, 3>> proteinCarbons;
        for (int i = 0; i < coordInfo.residues.size(); i++) {
            if (std::find(topOrder.begin(), topOrder.end(), coordInfo.residues[i]->res_name) == topOrder.end()) {
				for (const auto& weakAtom : coordInfo.residues[i]->atoms) {
					auto atom = weakAtom.lock();
					proteinAtomIndices.push_back(atom->atom_num + 1);
					if (atom->element == "C") {
						proteinCarbons.push_back(atom->coord);
					}
				}
            }
            else {
				break;
            }
        }

		// Get indices of non-protein, non-atmosphere residues 
		std::unordered_map<std::string, std::vector<int>> dropletIndices;
        std::unordered_map<std::string, std::vector<int>> atmosphereIndices;
        for (const auto& resType : topOrder) {
            dropletIndices[resType] = {};
            atmosphereIndices[resType] = {};

			if (coordInfo.residueMap[resType].empty()) {
				continue;
			}

            // N2 and O2 always gas
            if (resType == "NNN" || resType == "OOO") {
                for (int i = 0; i < coordInfo.residueMap[resType].size(); i++) {
                    atmosphereIndices[resType].push_back(i);
                }
            }

			// For others, MUST differentiate between droplet and atmospheric molecules!
            else {
                float cutoff = 5.0f; // nm - TODO - this may need to be increased for larger droplets
                for (int i = 0; i < coordInfo.residueMap[resType].size(); i++) {
                    auto firstAtom = coordInfo.residueMap[resType][i]->atoms[0].lock();
                    std::vector<float> dist = norm(proteinCarbons, firstAtom->coord);
                    if (*std::min_element(dist.begin(), dist.end()) > cutoff) {
                        atmosphereIndices[resType].push_back(i);
					}
                    else {
                        dropletIndices[resType].push_back(i);
                    }
                }
            }
        }

		// Begin writing .ndx file
        std::ofstream ndx_file(ndx_fname + ".ndx");
        if (!ndx_file.is_open()) {
			std::cerr << "Failed to open file: " << ndx_fname << ".ndx" << std::endl;
            std::exit(1);
        }

        // Write protein grp first
        ndx_file << "[ Protein ]\n";
        int count = 1;
        for (const auto& index : proteinAtomIndices) {
            ndx_file << std::right << std::setw(4) << index << " ";
            count++;
            if (count % 16 == 0) {
                ndx_file << '\n';
                count = 1;
            }
        }
		if (count != 1) {
			ndx_file << '\n';
		}

		// Finish by writing droplet_protein, droplet and atmosphere grps, record whether groups are actually present
		bool foo = writeNDXGroup(ndx_file, topOrder, "Protein_Droplet", dropletIndices, coordInfo.residueMap, true, proteinAtomIndices);
        bool droplet = writeNDXGroup(ndx_file, topOrder, "Droplet", dropletIndices, coordInfo.residueMap, false, {});
        bool atmosphere = writeNDXGroup(ndx_file, topOrder, "Atmosphere", atmosphereIndices, coordInfo.residueMap, false, {});
        ndx_file.close();

		// Now begin modiyfing .mdp file, start with comm_grps
		std::string commGrps = "Protein";
		if (droplet) {
			commGrps += "_Droplet";
		}
		if (atmosphere) {
			commGrps += " Atmosphere ";
		}

		// tc_grps
        std::string tcGrps = "Protein ";
        std::string refT = std::to_string(temperature) + " ";
        std::string tauT = "100 ";
		if (droplet) {
			tcGrps += "Droplet ";
			refT += std::to_string(temperature) + " ";
			tauT += "5 ";
		}
        if (atmosphere) {
            tcGrps += "Atmosphere ";
            refT += std::to_string(config.gas_temp) + " ";
            tauT += "5 ";
        }

        // Determine proper fmm-override-tree-depth if using FMM
        // TODO - this needs to be optimized
        std::string fmmDepth = "0";
        int numAtoms = coordInfo.atoms.size();
        if (numAtoms > 1000000) {
            fmmDepth = "4";
        }
        else if (numAtoms > 30000) {
            fmmDepth = "3";
        }
        else if (numAtoms > 20000) {
            fmmDepth = "2";
        }
        else if (numAtoms > 10000) {
            fmmDepth = "1";
        }
        else {
            fmmDepth = "0";
        }

        // Modify non-bonded cutoffs to fit within box
		std::string newCutoff = std::to_string((coordInfo.box_vectors[0]*0.99f) / 3.0f);

        // Modify .mdp file
        std::ifstream mdp_file(mdp_fname + ".mdp");
        if (!mdp_file.is_open()) {
            std::cerr << "Failed to open file: " << mdp_fname << ".mdp" << std::endl;
            return;
        }
        std::string line;
        std::vector<std::string> lines;
        while (std::getline(mdp_file, line)) {
            lines.push_back(line);
        }

        for (size_t i = 0; i < lines.size() - 1; ++i) {
            if (lines[i].find("comm_grps") != std::string::npos) {
                lines[i] = "comm_grps          = " + commGrps;
            }
            else if (lines[i].find("tc_grps") != std::string::npos) {
                lines[i] = "tc_grps            = " + tcGrps;
            }
            else if (lines[i].find("ref_t") != std::string::npos) {
                lines[i] = "ref_t              = " + refT;
            }
            else if (lines[i].find("tau_t") != std::string::npos) {
                lines[i] = "tau_t              = " + tauT;
            }
            else if (lines[i].find("fmm-override-tree-depth") != std::string::npos) {
                lines[i] = "fmm-override-tree-depth = " + fmmDepth;
            }
            else if (lines[i].find("rlist") != std::string::npos) {
                lines[i] = "rlist              = " + newCutoff;
            }
            else if (lines[i].find("rcoulomb") != std::string::npos) {
                lines[i] = "rcoulomb           = " + newCutoff;
            }
            else if (lines[i].find("rvdw") != std::string::npos) {
                lines[i] = "rvdw               = " + newCutoff;
            }
        }

        std::ofstream out(mdp_fname + ".mdp");
        for (const auto& line : lines) {
            out << line << '\n';
        }
    }

}