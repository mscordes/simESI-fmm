#include "Core.h"
#include <cstring>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <unordered_map>

namespace Core {

	// Master run loop for the simulation
	void simulation(const Config& config, const std::filesystem::path& trialPath) {

		// Order of residues in the top file, order must always be maintained so referenced frequently
		const std::vector<std::string> topOrder = {
			"SOL", "HHO", "OHX", "ATX", "AHX", "NXX", "NXH", "NNN", "OOO"
		};

		// Begin removing contaminants from .pdb, begin by creating .ndx file
		{
			std::ostringstream command;
			command << config.gmx_env << " make_ndx -f " << config.pdb << " -o clean.ndx";
			std::vector<std::string> inputs = { "q" };
			auto_gmx_input(command.str(), inputs);
		}

		// Use .ndx via editconf to actually remove non-Protein atoms
		{
			std::ostringstream command;
			command << config.gmx_env << " editconf -f " << config.pdb << " -n clean.ndx -o clean.pdb -box " <<
				std::to_string(config.box_size) << " -c";
			std::vector<std::string> inputs = { "1", "1" };
			auto_gmx_input(command.str(), inputs);
		}

		// Unpack clean .pdb file
		CoordInfo coordInfo = buildCoordInfo("clean.pdb", config.box_size);

		// Ensure correct ordering of atoms
		reorderAtoms(coordInfo, topOrder);

		// Write initial .pdb file with corrected residue/protein ordering
		writePDB("init.pdb", coordInfo.proteinAtoms, coordInfo.box_vectors);	
		coordInfo = buildCoordInfo("init.pdb", config.box_size);

		// Get pKa values using corrected .pdb
		std::unordered_map<int, float> pkaMap = getPkaVals(coordInfo, "init");

		// Droplet formation and equilibriation, check if continuing
		int initStep;
		if (config.dir_cont == "no") {
			Core::equilibriation(config, coordInfo, pkaMap, config.box_size, topOrder);
			copyFile(std::filesystem::path ("nvt.gro"), std::filesystem::path("0.gro"));
			copyFile(std::filesystem::path("droplet.top"), std::filesystem::path("0.top"));
			coordInfo = buildCoordInfo("0.gro", config.box_size);
			initStep = 1;
		}
		// If continuing from previous run
		else {
			if (config.step_cont == -1) {
				std::cerr << "Continuation specified without -step_cont argument." << std::endl;
				exit(1);
			}
			std::cout << "Continuing from previous simulation." << std::endl;
			initStep = config.step_cont;
			int failStep =  initStep - 1; 
			std::string gro = std::to_string(failStep) + ".gro";
			std::string top = std::to_string(failStep) + ".top";
			topFromCoord(gro, top, config.box_size, config, topOrder);
			coordInfo = buildCoordInfo(gro, config.box_size);
		}

		// Set water_cutoff arg based on protein mass if not user defined
		int water_cutoff = config.water_cutoff;
		if (config.water_cutoff == -1) {
			float proteinMass = getProteinMass(coordInfo.proteinAtoms);
			water_cutoff = static_cast<int>(std::round(0.0035f * proteinMass));
		}

		// Get initial net protein and droplet charge for writing to outputs
		std::vector<float> protCharges = getProtCharges(coordInfo.proteinAtoms);
		int protCharge = getNetCharge(protCharges);
		std::vector<float> charges = getCharges(coordInfo.proteinAtoms, coordInfo.numResidues, topOrder);
		int sysCharge = getNetCharge(charges);

		// Write initial outputs
		std::filesystem::path dataPath = trialPath / "data";
		for (const auto& resType : topOrder) {
			if (coordInfo.numResidues.find(resType) != coordInfo.numResidues.end()) {
				write_output(dataPath, std::to_string(coordInfo.numResidues.at(resType)), resType + ".txt");
			}
		}
		write_output(dataPath, std::to_string(protCharge), "prot_charge.txt");
		write_output(dataPath, std::to_string(sysCharge), "sys_charge.txt");

		// For GROMACS checkpoint file naming
		int numRestarts = 1;
		int numContSteps = 2000;

		// Initialize timers
		std::chrono::time_point<std::chrono::high_resolution_clock> initTime = timer();
		std::chrono::time_point<std::chrono::high_resolution_clock> initStepTime;
		std::chrono::time_point<std::chrono::high_resolution_clock> endStepTime;
		std::chrono::time_point<std::chrono::high_resolution_clock> gromppTime;
		std::chrono::time_point<std::chrono::high_resolution_clock> mdTime;
		std::chrono::time_point<std::chrono::high_resolution_clock> endTime;

		// Simulation loop
		int numSteps = static_cast<int>(std::round(config.time * 250)); // 1 ns = 250 4ps steps
		float temperature; // Droplet temperature
		for (int step = initStep; step < numSteps; step++) {
			initStepTime = timer();

			// Copy over coord/top files from previous step and build CoordInfo
			std::string oldStep = std::to_string(step - 1);
			std::string oldTop = oldStep + ".top";
			std::string oldGro = oldStep + ".gro";
			std::string newTop = std::to_string(step) + ".top";
			std::string newGro = std::to_string(step) + "_preMD.gro";
			copyFile(std::filesystem::path(oldTop), std::filesystem::path(newTop));
			coordInfo = buildCoordInfo(oldGro, config.box_size);

			// Recalculate pKa's to account for conformationl changes during droplet evaporation
			if (step % 250 == 0) {
				updatePkavals(pkaMap, coordInfo);
			}

			// Set droplet temperature based on user defined args
			temperature = (coordInfo.numResidues["SOL"] > config.water_cutoff) ? config.init_temp : config.final_temp;

			// Sometimes pdb2gmx will miss disulfide bonds if bond slightly too long/short, so correct here
			coordInfo.atoms = fixDisulfides(coordInfo.residueMap, coordInfo.atoms);

			// Flags to determine if need new protein top (prot) or only needs updating coordInfo (pairs)
			bool pairs = false;
			bool prot = false;
			bool createRunFile = false;

			/* Begin finding exchanges, repeated 5x for Grotthuss diffusion see,  
			Lars Konermann and Scott Kim Journal of Chemical Theory and Computation 2022 18 (6), 3781 - 3794 */
			std::vector<int> clusterIDs;
			std::unordered_map<int, Cluster> clusters;
			TitratableSites titSites;
			std::vector<std::array <float, 3>> skipCoords; // Prevent 'rattling' of proton between donor/acceptor pair
			for (int hop = 0; hop < 5; hop++) {

				// If no exchanges occurred past first hop, break
				if (hop > 0 && !pairs) {
					break;
				}

				// Find cluster and titratable site information on first hop, or if exchanges have occurred
				if (hop == 0 || pairs) {

					// Find the cluster each water belongs to via Kruskal's algorithm
					clusterIDs = computeClusters(coordInfo.waterOCoords);

					// Get number of waters, ions, and pH of each cluster
					clusters = binClusters(clusterIDs, coordInfo, config);

					// Get titratable site information
					titSites = getTitratableSites(coordInfo, clusterIDs, clusters, pkaMap);

					// Get partial charges
					charges = getCharges(coordInfo.proteinAtoms, coordInfo.numResidues, topOrder);
				}

				// Reset flags
				pairs = false;
				prot = false;

				// Find proton donor/acceptor pairs, if restarting avoid exchanges for 1 step
				std::vector<Exchange> exchanges = findPairs(coordInfo, titSites, clusterIDs,
					clusters, temperature, step, hop, charges, skipCoords, pairs, prot);

				// If continuing from failure, skip exchanges for 1 step
				if (config.dir_cont != "no" && step == config.step_cont) {
					exchanges.clear();
				}

				// Write accepted exchanges to output 
				for (const auto& exchange : exchanges) {
					write_output(dataPath, exchange.toString(), "exchanges.txt");
				}

				// Facilitate proton exchanges, and update coordInfo if necessary
				doExchanges(coordInfo, exchanges, pairs, prot, topOrder, config);

				// Update coord/top files, if protein top needs updating, write new .top
				if (prot) {
					std::unordered_map<std::string, std::vector<std::string>> pdb2gmxMap = getProtStates(coordInfo);
					writeTOP(newTop, coordInfo, pdb2gmxMap, false, config, topOrder);
					createRunFile = true;
				}
				// Else if only non-protein exchanges, just update the tally of each residue in the .top
				else if (pairs) {
					createTOP(newTop, coordInfo, topOrder);
					createRunFile = true;
				}
			}

			// Remove evaporated atoms
			std::vector<std::array<float, 3>> proteinCarbons;
			std::unordered_map<std::string, int> vaporDict;
			bool evap = removeEvaporated(config, coordInfo, topOrder, newTop, proteinCarbons, vaporDict);

			// Recenter droplet to prevent interactions with PBCs
			bool recentered = recenterDroplet(coordInfo, proteinCarbons);

			// Gas velocities can become distorted due to stitching, correct here
			coordInfo.residues = binAtoms(coordInfo.atoms);
			coordInfo.residueMap = getResidueMap(coordInfo.residues);
			bool gasCorrect = correctGasVelocities(coordInfo.atoms, temperature);
			fixAmmoniaGas(config, coordInfo, gasCorrect); // Fix instability bug with gaseous ammonia
			fixAceticGas(config, coordInfo, gasCorrect); // Fix instability bug with gaseous acetic acid

			endStepTime = timer();

			// If molecules have changed, create new run file and run MD
			if (createRunFile || evap || recentered || gasCorrect || step == 1 || step % 50 == 0 || step == config.step_cont ) {

				// Write final post-exchange .gro file
				writeGRO(newGro, coordInfo.atoms, coordInfo.box_vectors);

				// Prep .ndx and .mdp files for MD, tc-grps will break unless the specified residue name is present
				modifyMDPgrps(std::to_string(step) + "_preMD", "prodrun", config, topOrder, coordInfo, temperature);

				// Run MD
				{
					std::ostringstream command;
					command << config.gmx_env << " grompp -f prodrun.mdp -c " << newGro << " -p " << newTop << " -n "
						<< std::to_string(step) + "_preMD.ndx -o " << std::to_string(step) << ".tpr -maxwarn 100";
					std::vector<std::string> inputs = {};
					auto_gmx_input(command.str(), inputs);
					gromppTime = timer();
				}
				deleteFile(std::to_string(step) + ".gro");
				runMD(std::to_string(step), config);
				mdTime = timer();

				// Update checkpoint naming vars
				numRestarts = 1;
				numContSteps = 2000;

				// Clean up outputs
				deleteFile(std::to_string(step) + ".edr");
				deleteFile(std::to_string(step) + ".log");

				createRunFile = false;
			}

			// If no molecules have changed, use GROMACS checkpointing
			else {

				// Update checkpoint naming vars
				numRestarts++;
				numContSteps += 2000;
				copyFile(std::filesystem::path(std::to_string(step - 1) + ".cpt"), std::filesystem::path(std::to_string(step) + ".cpt"));

				// Extend time in run file
				{
					std::ostringstream command;
					command << config.gmx_env << " convert-tpr -s " << std::to_string(step - 1) << ".tpr -o " 
						<< std::to_string(step) << ".tpr -nsteps " << std::to_string(numContSteps);
					std::vector<std::string> inputs = {};
					auto_gmx_input(command.str(), inputs);
					gromppTime = timer();
				}

				// Run MD
				{
					std::ostringstream command;
					command << config.gmx_env << " mdrun -cpi " << std::to_string(step) << " -deffnm " 
						<< std::to_string(step) << " -noappend -cpo " << std::to_string(step);
					// Modify based on resources
					if (config.hpc == "yes") {
						command << " -ntmpi 1 -nb gpu";
					}
					else if (config.gpu == "yes") {
						command << " -nb gpu";
					}
					std::vector<std::string> inputs = {};
					auto_gmx_input(command.str(), inputs);
					mdTime = timer();
				}

				// Rename outputs to account for GROMACS checkpointing
				for (const auto& ftype : { ".gro", ".xtc", ".edr", ".log" }) {
					renameCheckpointFile(ftype, step, numRestarts);
				}
			}
			 
			// Check if MD completed successfully
			if (!std::filesystem::exists(std::to_string(step) + ".gro")) {
				std::cerr << "MD failed to complete successfully." << std::endl;
				exit(1);
			}

			// Clean up files 
			deleteFile(std::filesystem::path(std::to_string(step) + ".edr"));
			deleteFile(std::filesystem::path(std::to_string(step) + ".log"));
			deleteFile(std::filesystem::path(newGro));

			// Save every 1 of every 10 traj files
			if (step % 10 != 0) {
				deleteFile(std::filesystem::path(std::to_string(step) + ".xtc"));
			}
			// Delete older other files, need to save for potential checkpointing
			if (step > 2) {
				int older_step = step - 2;
				for (const auto& ftype : { ".cpt", ".tpr", ".top" }) {
					deleteFile(std::filesystem::path(std::to_string(older_step) + ftype));
				}
				deleteFile(std::filesystem::path(std::to_string(older_step) + "_prev.cpt"));
				deleteFile(std::filesystem::path(std::to_string(older_step) + "_preMD.ndx"));
			}

			// Write outputs, start with how much of each solute
			for (const auto& resType : topOrder) {
				if (coordInfo.numResidues.find(resType) != coordInfo.numResidues.end()) {
					write_output(dataPath, std::to_string(coordInfo.numResidues.at(resType)), resType + ".txt");
				}
			}

			// Write how much each gas phase solute present
			for (const auto& [resType, num] : vaporDict) {
				if (resType != "NNN" && resType != "OOO") {
					write_output(dataPath, std::to_string(num), "vapor_" + resType + ".txt");
				}
			}

			// System charge and temperature
			write_output(dataPath, std::to_string(getNetCharge(charges)), "sys_charge.txt");
			protCharges = getProtCharges(coordInfo.proteinAtoms);
			write_output(dataPath, std::to_string(getNetCharge(protCharges)), "prot_charge.txt");
			write_output(dataPath, std::to_string(temperature), "temperature.txt");

			// Timing information
			auto endTime = timer();
			write_output(dataPath, duration(initStepTime, endStepTime), "exchange_time.txt");
			write_output(dataPath, duration(endStepTime , gromppTime), "grompp_time.txt");
			write_output(dataPath, duration(gromppTime  , mdTime), "md_time.txt");
			write_output(dataPath, duration(initTime    , endTime), "tot_time.txt");

			std::cout << "\n\n\n########################################" << std::endl;
			std::cout << "          STEP " << step << " COMPLETED          " << std::endl;
			std::cout << "########################################\n\n\n" << std::endl;
		}
		std::cout << "Time limit of " << config.time << " ns reached." << std::endl;
		std::cout << "Simulation complete, final charge is " << getNetCharge(charges) << "." << std::endl;
		std::exit(0);
	}
}