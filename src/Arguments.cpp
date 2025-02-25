#include "Core.h"
#include <iostream>
#include <string>
#include <stdexcept>
#include <unordered_map>

namespace Core {

    void printHelp() {
        std::cout << "Usage: ./parser -pdb <filename> [options]\n\n"
            << "Required arguments:\n"
            << "  -pdb <filename>        PDB filename in simESI-main/input_files.\n\n"
            << "Optional arguments:\n"
            << "  -esi_mode <pos|neg>    Positive ('pos') or negative ('neg') ion mode. Default: 'pos'.\n"
			<< "  -atm <yes|no>          Simulate with or without atmosphere. Default of \"yes\".\n"
            << "  -amace_conc <float>    Initial ammonium acetate concentration. Default: 0.25.\n"
            << "  -water_vapor <float>   Percent mass of water vapor to seed in atmosphere. Default: 0 (dry atmosphere).\n"
            << "  -ach_vapor <float>     Percent mass of acetic acid vapor to seed in atmosphere. Default: 0.\n"
            << "  -nh3_vapor <float>     Percent mass of ammonia vapor to seed in atmosphere. Default: 0.\n"
            << "  -ace_vapor <float>     (EXPERIMENTAL) Percent mass of acetate vapor to seed in atmosphere. Default: 0.\n"
            << "  -nh4_vapor <float>     (EXPERIMENTAL) Percent mass of ammonium vapor to seed in atmosphere. Default: 0.\n"
            << "  -droplet_size <float>  How much water to seed around protein in nm, default of 1.5 nm.\n"
            << "  -box_size <float>      Simulation box size, ie, how much atmosphere. Default of 100 nm (^3).\n"
            << "  -time <float>          Simulation cutoff time in ns. Default: 25.\n"
            << "  -init_temp <float>     Initial droplet temperature in K. Default: 370.\n"
            << "  -final_temp <float>    Final droplet temperature to completely desovlate. Default: 450.\n"
            << "  -gas_temp <float>      Temperature of the surrounding atmosphere in K. Default: 300.\n"
            << "  -water_cutoff <int>    Remaining number of waters to ramp temperature." << 
                                            "Default calculated based on protein mass\n"
            << "  -pka_pdb <string>      Alternative .pdb for PROPKA pka value determination.\n"
            << "  -dir_cont <dir>        Continue run, input dirname of previous trial to continue from.\n"
            << "  -step_cont <int>       If continuing, step to continue from.\n"
            << "  -gmx_env <string>      Path to gmx executable, default is \"gmx\" assuming gmx is in path.\n"
            << "  -gpu <yes|no>          Use GPU acceleration. Default \"yes\", CPU only with \"no\".\n"
            << "  -hpc <yes|no>          Modify gmx calls for HPC use. Default \"no\".\n"
            << "  -fmm <yes|no>          Use FMM enabled gmx. Default \"no\".\n";
    }

    void validateConfig(const Config& config) {
        if (config.pdb.empty()) {
            std::cerr << "Missing required argument: -pdb.";
            exit(1);
        }
        if (config.esi_mode != "pos" && config.esi_mode != "neg") {
            std::cerr << "Invalid value for -esi_mode. Must be 'pos' or 'neg'.";
            exit(1);
        }
        if (config.atm != "yes" && config.atm != "no") {
            std::cerr << "Invalid value for -atm. Must be 'yes' or 'no'.";
            exit(1);
        }
        if (config.gpu != "yes" && config.gpu != "no") {
            std::cerr << "Invalid value for -gpu. Must be 'yes' or 'no'.";
            exit(1);
        }
        if (config.hpc != "yes" && config.hpc != "no") {
            std::cerr << "Invalid value for -hpc. Must be 'yes' or 'no'.";
            exit(1);
        }
        if (config.fmm != "yes" && config.fmm != "no") {
            std::cerr << "Invalid value for -fmm. Must be 'yes' or 'no'.";
            exit(1);
        }
    }

    // Parse command line arguments
    Config getArgs(int argc, char* argv[]) {
        if (argc == 1 || std::string(argv[1]) == "-help") {
            printHelp();
            exit(0);
        }

        Config config;
        std::unordered_map<std::string, std::string> args;

        for (int i = 1; i < argc; i += 2) {
            if (i + 1 >= argc || argv[i][0] != '-') {
                std::cerr << "Invalid arguments. Use -help for usage.";
                exit(1);
            }
            args[argv[i]] = argv[i + 1];
        }

        // Extract and assign values
        for (const auto& arg : args) {
            if (arg.first == "-pdb") config.pdb = args["-pdb"];
            else if (arg.first == "-esi_mode") config.esi_mode = args["-esi_mode"];
            else if (arg.first == "-atm") config.atm = args["-atm"];
			else if (arg.first == "-amace_conc") config.amace_conc = std::stof(args["-amace_conc"]);
			else if (arg.first == "-water_vapor") config.water_vapor = std::stof(args["-water_vapor"]);
            else if (arg.first == "-ace_vapor") config.ace_vapor = std::stof(args["-ace_vapor"]);
            else if (arg.first == "-ach_vapor") config.ach_vapor = std::stof(args["-ach_vapor"]);
            else if (arg.first == "-nh4_vapor") config.nh4_vapor = std::stof(args["-nh4_vapor"]);
            else if (arg.first == "-nh3_vapor") config.nh3_vapor = std::stof(args["-nh3_vapor"]);
            else if (arg.first == "-droplet_size") config.droplet_size = std::stof(args["-droplet_size"]);
            else if (arg.first == "-box_size") config.droplet_size = std::stof(args["-box_size"]);
            else if (arg.first == "-time") config.time = std::stof(args["-time"]);
            else if (arg.first == "-init_temp") config.init_temp = std::stof(args["-init_temp"]);
            else if (arg.first == "-final_temp") config.final_temp = std::stof(args["-final_temp"]);
            else if (arg.first == "-gas_temp") config.gas_temp = std::stof(args["-gas_temp"]);
            else if (arg.first == "-water_cutoff") config.water_cutoff = std::stoi(args["-water_cutoff"]);
            else if (arg.first == "-pka_pdb") config.pka_pdb = args["-pka_pdb"];
            else if (arg.first == "-dir_cont") config.dir_cont = args["-dir_cont"];
            else if (arg.first == "-step_cont") config.step_cont = std::stoi(args["-step_cont"]);
            else if (arg.first == "-gmx_env") config.gmx_env = args["-gmx_env"];
            else if (arg.first == "-gpu") config.gpu = args["-gpu"];
            else if (arg.first == "-hpc") config.hpc = args["-hpc"];
            else if (arg.first == "-fmm") config.fmm = args["-fmm"];
            else {
				std::cerr << "Invalid argument: " << arg.first << std::endl; 
                exit(1);
            }
        }

        validateConfig(config);

        // Display parsed arguments
        std::cout << "PDB file: " << config.pdb << std::endl;
        std::cout << "ESI mode: " << config.esi_mode << std::endl;
        std::cout << "Under atmosphere: " << config.atm << std::endl;
        std::cout << "Ammonium acetate concentration: " << config.amace_conc << " M" << std::endl;
        std::cout << "Time cutoff: " << config.time << " ns" << std::endl;
        std::cout << "Initial temperature: " << config.init_temp << " K" << std::endl;
        std::cout << "Final temperature: " << config.final_temp << " K" << std::endl;
        std::cout << "Gas temperature: " << config.gas_temp << " K" << std::endl;
		std::cout << "gmx executable: " << config.gmx_env << std::endl;
		std::cout << "GPU acceleration: " << config.gpu << std::endl;
		std::cout << "HPC mode: " << config.hpc << std::endl;
        std::cout << "FMM enabled: " << config.fmm << std::endl;

        return config;
    }
}