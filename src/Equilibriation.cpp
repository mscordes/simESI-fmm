#include "Core.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <sys/stat.h>

namespace Core {

    // Sanitizes .pdb, creates droplet, and equilibriates
    void equilibriation(const Config& config, CoordInfo& coordInfo, const std::unordered_map<int, float>& pkaMap, const float& init_boxSize, 
        const std::vector<std::string>& topOrder) {

        /* Test validity of init.pdb and convert to a proper CHARMM36 formatted.gro file vid pdb2gmx
           First need to set inputs for termini for pdb2gmx */
        std::vector<std::string> inputs;
        inputs.reserve(coordInfo.proteins.size() * 2); // Reserve space to avoid reallocations
        for (const auto& monomer : coordInfo.proteins) {
            // N-termini
            auto nterm = monomer->residues[0].lock();
            inputs.push_back((nterm->res_name == "PRO" || nterm->res_name == "MET") ? "1" : "0");

            // C-termini always 0 in CHARMM36
            inputs.push_back("0");
        }

        // Call pdb2gmx
        {
            std::ostringstream command;
            command << config.gmx_env << " pdb2gmx -f init.pdb -o init.gro -p temp.top -ff charmm36 -water none -ignh -ter";
            auto_gmx_input(command.str(), inputs);
        }

        // Check if topology generation failed, if success, delete posre.itp file 
        {
            std::ifstream file("temp.top");
            if (!file.is_open()) {
                std::cerr << "\n\n\nERROR: Could not form molecular topology. This is likely due to not having"
                    << " the GROMACS executable 'gmx' callable from command line." << std::endl;
                std::exit(1);
            }
            std::string line;
            int line_count = 0;
            while (std::getline(file, line)) {
                if (++line_count > 1) {
                    break;
                }
            }
            if (line_count < 2) {
                std::cerr << "\n\n\nInputted pdb file is compromised in some way, unable to generate topology." << std::endl;
                std::exit(1);
            }

			// Delete posre.itp file
            deleteFile(std::filesystem::path("posre.itp"));
        }

        // Unpack CHARMM formatted .gro
        coordInfo = buildCoordInfo("init.gro", init_boxSize);

        // Get expected protonation states at pH 7
        const std::unordered_map<std::string, std::vector<std::string>> pdb2gmxMap = setProtStates(coordInfo, pkaMap, config);

        // Make "system.gro" and "system.top" files determined protonation states
        writeTOP("system.top", coordInfo, pdb2gmxMap, true, config, topOrder);

        // Update coordInfo 
        coordInfo = buildCoordInfo("system.gro", init_boxSize);

        // Reorder atoms to maintain topology
        reorderAtoms(coordInfo, topOrder);

        // Droplet formation
        coordInfo = formDroplet(config, coordInfo, topOrder);

        // Energy minimization
        {
            std::ostringstream command;
            command << config.gmx_env << " grompp -f em.mdp -c droplet.gro -p droplet.top -o em.tpr -maxwarn 100";
            std::vector<std::string> inputs = {};
            auto_gmx_input(command.str(), inputs);
            runMD("em", config);
        }

        // Modify tc-grps and comm-grps in nvt.mdp
        modifyNDXGrps("droplet", coordInfo.numResidues, config);
        modifyMDPGrps("nvt", coordInfo.numResidues, config.init_temp, coordInfo.atoms);

        // NVT equilibriation
        {
            std::ostringstream command;
            command << config.gmx_env << " grompp -f nvt.mdp -c em.gro -p droplet.top -n droplet.ndx -o nvt.tpr -maxwarn 100";
            std::vector<std::string> inputs = {};
            auto_gmx_input(command.str(), inputs);
            runMD("nvt", config);
        }

		// Check if equilibriation completed successfully
		{
			std::ifstream file("nvt.gro");
			if (!file.is_open()) {
				std::cerr << "\n\n\nERROR: NVT equilibriation failed." << std::endl;
				std::exit(1);
			}
            else {
                std::cout << "\n\n\nNVT equilibriation successful." << std::endl;
            }
		}
    }
}