#include "Core.h"
#include <sstream>
#include <fstream>
#include <iostream>

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

	// Create .ndx file for COM translational/rotational removal and T-coupling
    void modifyNDXGrps(const std::string& fname, const std::unordered_map<std::string, int>& numResidues, 
        const Config& config) {

		std::vector<std::string> inputs;
        std::string commGrps = "1 |";
		std::string tcGrps = "";

		std::vector<std::string> resTypes = { "SOL", "HHO", "OHX", "ATX", "AHX", "NXX", "NXH" };
        for (const auto& resType : resTypes) {
            auto it = numResidues.find(resType);
            if (it != numResidues.end() && it->second > 0) {
                commGrps += " r " + resType + " |";
                tcGrps += " r " + resType + " |";
            }
        }
		commGrps.pop_back();
        inputs.push_back(commGrps);

		if (!tcGrps.empty()) {
			tcGrps.pop_back();
			inputs.push_back(tcGrps);
		}

        inputs.push_back("r NNN | r OOO");
		inputs.push_back("q");


        std::ostringstream command;
        command << config.gmx_env << " make_ndx -f " << fname << ".gro -o " << fname << ".ndx";
        auto_gmx_input(command.str(), inputs);
    }

    // Modify .mdp file for COM translational/rotational removal and T-coupling
    void modifyMDPGrps(const std::string& fname, const std::unordered_map<std::string, int>& numResidues, const float& temperature,
        const std::vector<std::shared_ptr<Atom>> atoms) {

        std::string commGrps = "Protein_";
        std::string tcGrps = "Protein ";
        std::string refT = std::to_string(temperature) + " ";
        std::string tauT = "100 ";

        std::string new_tcGrps = "";
        std::vector<std::string> resTypes = { "SOL", "HHO", "OHX", "ATX", "AHX", "NXX", "NXH" };
        for (const auto& resType : resTypes) {
            auto it = numResidues.find(resType);
            if (it != numResidues.end() && it->second > 0) {
                commGrps += resType + "_";
                new_tcGrps += resType + "_";
            }
        }
        commGrps.pop_back();

        if (!new_tcGrps.empty()) {
            new_tcGrps.pop_back();
			tcGrps += new_tcGrps;
            refT += std::to_string(temperature) + " ";
            tauT += "5 ";
        }

		commGrps += " NNN_OOO ";
		tcGrps += " NNN_OOO ";
        refT += "300 ";
		tauT += "5 ";

        // Determine proper fmm-override-tree-depth
		std::string fmmDepth = "0";
        int numAtoms = atoms.size();
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

        // Modify .mdp file
		std::ifstream file(fname + ".mdp");
		if (!file.is_open()) {
			std::cerr << "Failed to open file: " << fname << ".mdp" << std::endl;
			return;
		}
        std::string line;
        std::vector<std::string> lines;
        while (std::getline(file, line)) {
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
        }

		std::ofstream out(fname + ".mdp");
		for (const auto& line : lines) {
			out << line << '\n';
		}
    }
}