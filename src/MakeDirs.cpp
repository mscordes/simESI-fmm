#include "Core.h"
#include <filesystem>

namespace Core {

	// Make and populate trial directory with necessary input files
	std::filesystem::path makeDirs(Config& config) {

		// Get simESI-cpp root
		std::filesystem::path currentPath = std::filesystem::current_path();
		std::filesystem::path rootPath = currentPath.parent_path();
		
		// Get coordinateFiles and inputFiles dirs
		std::filesystem::path coordinateFilesPath = rootPath / "coordinateFiles";
		std::filesystem::path inputFilesPath = rootPath / "inputFiles";
		if (!std::filesystem::exists(coordinateFilesPath) || !std::filesystem::exists(inputFilesPath)) {
			rootPath = currentPath;
			coordinateFilesPath = rootPath / "coordinateFiles";
			inputFilesPath = rootPath / "inputFiles";
		
			if (!std::filesystem::exists(coordinateFilesPath) || !std::filesystem::exists(inputFilesPath)) {
				std::cerr << "coordinateFiles or inputFiles directory does not exist in current directory." << std::endl;
				std::exit(1);
			}
		}
		std::cout << "simESI-root directory: " << rootPath << std::endl;
		
		// Make outputFiles dir if not exists
		std::filesystem::path outputFilesPath = rootPath / "outputFiles";
		if (!std::filesystem::exists(outputFilesPath)) {
			std::filesystem::create_directory(outputFilesPath);
		}
		
		// Remove .pdb from the inputted .pdb path
		std::string pdbName;
		if (config.pdb.length() >= 5) {
			// Use substr to get the substring from the beginning to the length minus 4
			pdbName = config.pdb.substr(0, config.pdb.length() - 4);
		}
		else {
			std::cout << "Inputted .pdb name too short, maybe forgot to end with \".pdb\"?" << std::endl;
			std::exit(1);
		}

		// Find/make trial dir as well as nested data and simulation dirs
		std::filesystem::path trialPath;
		std::filesystem::path dataPath;
		std::filesystem::path simulationPath;
		
		// Make if not continuation
		if (config.dir_cont == "no") {
			// Count up trial number 
			int trial_num = 1;
			while (std::filesystem::exists(outputFilesPath / (pdbName + "_" + std::to_string(trial_num)))) {
				trial_num++;
			}
			trialPath = outputFilesPath / (pdbName + "_" + std::to_string(trial_num));

			std::filesystem::create_directory(trialPath);
			dataPath = trialPath / "data";
			std::filesystem::create_directory(dataPath);
			simulationPath = trialPath / "simulation";
			std::filesystem::create_directory(simulationPath);
		}
		// Check if exists if continuation
		else {
			trialPath = outputFilesPath / config.dir_cont;
			if (!std::filesystem::exists(trialPath)) {
				std::cerr << "Continuation directory does not exist." << std::endl;
				std::exit(1);
			}
			dataPath = trialPath / "data";
			simulationPath = trialPath / "simulation";
		}

		// Begin populating simulation dir, start by copying .pdb
		std::filesystem::path pdbPath = coordinateFilesPath / config.pdb;
		if (!std::filesystem::exists(simulationPath / config.pdb) || std::filesystem::exists(pdbPath)) {
			copyFile(pdbPath, simulationPath / config.pdb);
		}
		else if (!std::filesystem::exists(pdbPath)) {
			std::cerr << "PDB file does not exist in coordinateFiles." << std::endl;
			std::exit(1);
		}

		// Copy forcefield dir
		std::filesystem::path ffPath = inputFilesPath / "charmm36.ff";
		if (!std::filesystem::exists(simulationPath / "charmm36.ff")) {
			copyRecursive(ffPath, simulationPath / "charmm36.ff");
		}

		// Copy necessary .mdp and .top files, selects .mdp's depending on if using FMM
		std::filesystem::path mdpPath = inputFilesPath / "mdp_files";
		copyFile(mdpPath / "em.mdp", simulationPath / "em.mdp");
		if (config.fmm == "yes") {
			copyFile(mdpPath / "nvt_FMM.mdp", simulationPath / "nvt.mdp");
			copyFile(mdpPath / "prodrun_FMM.mdp", simulationPath / "prodrun.mdp");
		}
		else {
			copyFile(mdpPath / "nvt.mdp", simulationPath / "nvt.mdp");
			copyFile(mdpPath / "prodrun.mdp", simulationPath / "prodrun.mdp");
		}
		std::filesystem::path topPath = inputFilesPath / "top_files";
		copyFiles(topPath, simulationPath);

		// Change to simulation dir
		try {
			std::filesystem::current_path(simulationPath);
			std::cout << "Changed to simulation dir: " << std::filesystem::current_path() << std::endl;
		}
		catch (const std::filesystem::filesystem_error& e) {
			std::cerr << "Error changing directory: " << e.what() << std::endl;
		}

		return trialPath;
	}
}