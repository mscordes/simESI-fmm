#include "Core.h"
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <random>
#include <set>
#include <cstring>
#include <algorithm>
#include <memory>
#include <filesystem>

namespace Core {

	float getProteinMass(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms) {

		// Atomic masses
		const std::unordered_map<std::string, float> atomMasses = {
			{"H", 1.008f}, {"C", 12.011f}, {"N", 14.007f}, {"O", 15.999f}, {"S", 32.06f}, {"M", 0.00f}
		};

		// Calculate protein mass
		float protMass = 0.0;
		for (const auto& [chain, atoms] : proteinAtoms) {
			for (const auto& atom : atoms) {
				float mass = atomMasses.at(atom->element);
				protMass += mass;
			}
		}

		return protMass;
	}

	// Carve out droplet around protein
	static std::vector<std::shared_ptr<Atom>> carveDroplet(const CoordInfo& coordInfo, const float& boxSize, const float& dropletRadius, 
		const std::vector<std::array<float, 3>>& proteinCoords) {

		std::vector<std::shared_ptr<Atom>> dropletAtoms;

		const std::array<float, 3> center = { boxSize / 2.0f, boxSize / 2.0f, boxSize / 2.0f };

		for (const auto& residue : coordInfo.residues) {
			if (residue->res_name == "SOL" && residue->atoms[0].lock()->atom_name == "OW") {
				std::vector<float> dist = norm(proteinCoords, residue->atoms[0].lock()->coord);
				if (*std::min_element(dist.begin(), dist.end()) < dropletRadius) {
					for (const auto& atom : residue->atoms) {
						dropletAtoms.push_back(atom.lock());
					}

				}
			}
			else {
				for (const auto& atom : residue->atoms) {
					dropletAtoms.push_back(atom.lock());
				}
			}
		}

		return dropletAtoms;
	}

	static std::array<float, 3> generateRandomLoc(const float& box_center, const float& droplet_radius) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> dis(static_cast<float>(box_center - (droplet_radius * 0.9f)),
			static_cast<float>(box_center + (droplet_radius * 0.9f)));

		// Generate random coordinates
		std::array<float, 3> random_loc = { dis(gen), dis(gen), dis(gen) };

		return random_loc;
	}

	// Insert a molecule of a given type into the droplet
	static void insertMolecule(const std::string& resGro, const int& numNew, std::vector<std::shared_ptr<Atom>>& dropletAtoms,
		std::vector<std::array<float, 3>>& coords, const float& dropletRadius, const std::array<float, 3>& center, 
		const float& minSeperation, const std::vector<std::array<float, 3>>& proteinCoords, const float& envelope) {

		// New molecule atoms and coords
		std::vector<std::shared_ptr<Atom>> molecAtoms = readGROatoms(resGro);
		std::vector<std::array<float, 3>> molecCoords = extractCoordinates(molecAtoms);

		// Begin inserting molecules randomly within droplet
		float minSepSquared = minSeperation * minSeperation;
		int resNum = dropletAtoms.back()->res_num;
		for (int i = 0; i < numNew; i++) {
			int attempts = 0; // Insertion attempts
			bool inserted = false;
			while (!inserted) {
				// Randomly select a location for the new molecule
				std::array<float, 3> newLoc = generateRandomLoc(center[0], dropletRadius);

				// Check is within enevlope around protein
				bool inDroplet = true;
				std::vector<float> dist = norm(proteinCoords, newLoc);
				if (*std::min_element(dist.begin(), dist.end()) > envelope) {
					inDroplet = false;
					continue;
				}

				// Molecule coordinates
				std::vector<std::array<float, 3>> newCoords;
				for (const auto& coord : molecCoords) {
					std::array<float, 3> newCoord = { coord[0] + newLoc[0], coord[1] + newLoc[1], coord[2] + newLoc[2] };
					newCoords.push_back(newCoord);
				}

				// Check if the molecule is too close to another molecule
				bool tooClose = false;
				for (const auto& newCoord : newCoords) {
					for (const auto& coord : coords) {
						const float dx = newCoord[0] - coord[0];
						const float dy = newCoord[1] - coord[1];
						const float dz = newCoord[2] - coord[2];
						const float distanceSquared = dx * dx + dy * dy + dz * dz;

						if (distanceSquared < minSepSquared) {
							tooClose = true;
							break;
						}
					}
				}


				// Insert the molecule into the droplet
				if (inDroplet && !tooClose) {
					resNum++;
					for (size_t i = 0; i < molecAtoms.size(); i++) {
						coords.push_back(newCoords[i]);

						std::shared_ptr<Atom> newAtom = std::make_shared<Atom>(*molecAtoms[i]);
						for (int j = 0; j < 3; j++) {
							newAtom->coord[j] += newLoc[j];
						}
						newAtom->res_id = newAtom->res_name + std::to_string(resNum);
						dropletAtoms.push_back(newAtom);
					}
					inserted = true;
				}

				// If too many attempts, break
				attempts++;
				if (attempts > 10000) {
					std::cerr << "ERROR: Could not insert all molecules during droplet creation." << std::endl;
					std::exit(1);
				}
			}
		}
	}

	// Complete droplet formation
	CoordInfo formDroplet(const Config& config, CoordInfo& coordInfo, const std::vector<std::string>& topOrder) {

		// Center protein in box just large enough to fit droplet
		std::array<float, 3> center = { 
			coordInfo.box_vectors[0] / 2.0f, 
			coordInfo.box_vectors[1] / 2.0f, 
			coordInfo.box_vectors[2] / 2.0f 
		};
		std::vector<std::array<float, 3>> proteinCoords = getProteinCoords(coordInfo.proteinAtoms);
		std::vector<float> dist = norm(proteinCoords, center);
		float newBoxSize = (*std::max_element(dist.begin(), dist.end()) + config.droplet_size) * 2.0f * 1.1f;
		{
		    std::ostringstream command;
		    command << config.gmx_env << " editconf -f system.gro -o system.gro -box " <<
		        std::to_string(newBoxSize) << " -c";
		    std::vector<std::string> inputs {};
		    auto_gmx_input(command.str(), inputs);
		}
		coordInfo = buildCoordInfo("system.gro", newBoxSize);

		// Find number of waters accounting for displacement by protein, start by solvating
		copyFile(std::filesystem::path("system.top"), std::filesystem::path("temp.top"));
		{
			std::ostringstream command;
			std::vector<std::string> inputs = {};
			command << config.gmx_env << " solvate -cp system.gro -cs tip4p_2005.gro -p temp.top -o temp.gro";
			auto_gmx_input(command.str(), inputs);
		}
		CoordInfo temp_coordInfo = buildCoordInfo("temp.gro", newBoxSize);

		// Carve droplet
		proteinCoords = getProteinCoords(temp_coordInfo.proteinAtoms);
		std::vector<std::shared_ptr<Atom>> tempAtoms = carveDroplet(temp_coordInfo, newBoxSize, config.droplet_size, proteinCoords);
		writeGRO("temp.gro", tempAtoms, temp_coordInfo.box_vectors);
		temp_coordInfo = buildCoordInfo("temp.gro", newBoxSize);
		int numSol = static_cast<int>(temp_coordInfo.residueMap.at("SOL").size());
		float solMass = numSol * 18.0f; // Water mass

		// Get radius of droplet assuming sphereical vol
		float protMass = getProteinMass(coordInfo.proteinAtoms);
		float protVol = protMass * 0.001361f; // Protein volume (nm^3)
		float solVol = solMass * 0.001660f; // Water volume (nm^3)
		float dropletVol = protVol + solVol;
		float dropletRadius = std::pow((3.0f / (4.0f * 3.14159f)) * dropletVol, 1.0f / 3.0f); // nm

		// Calculate 90% of Rayleigh limit as initial droplet net charge
		float rayleighConst = 1.258677f * std::pow(10.0f, 14.0f); // (8 * pi / e) * (eo * gamma)^(1/2)
		float rayleighRadius = std::pow(dropletRadius * std::pow(10.0f, -9.0f), 3.0f / 2.0f); // radius in m^(3/2)
		float rayleighLim = rayleighConst * rayleighRadius;
		float dropletCharge = static_cast<int>(std::round(0.9f * rayleighLim));

		// Get charge of solute to seed accounting for protein charge
		std::vector<float> protCharges = getProtCharges(coordInfo.proteinAtoms);
		int protCharge = getNetCharge(protCharges);
		int soluteCharge = dropletCharge - std::abs(protCharge);

		// Begin finding how many of each residue type to add now that we have number of waters
		std::unordered_map<std::string, int> dropletResidues{};

		// Number of base ammonium and acetate to seed (depending on inputted AmAc concentration)
		float numAmAc = 0.01801f * numSol * config.amace_conc;

		// Excess AmAc dependent on ESI mode (for net droplet charge), also effect from pH
		if (config.esi_mode == "pos") {
			dropletResidues["ATX"] = static_cast<int>(std::round(numAmAc*0.9));				// Acetate
			dropletResidues["AHX"] = static_cast<int>(std::round(numAmAc * 0.1));			// Acetic acid	
			int correction = static_cast<int>(std::round(numAmAc)) - static_cast<int>(std::round(numAmAc * 0.9)); // Correction HH.
			dropletResidues["NXH"] = soluteCharge + static_cast<int>(std::round(numAmAc)) - correction;	// Ammonium
		}
		else if (config.esi_mode == "neg") {
			int correction = static_cast<int>(std::round(numAmAc)) - static_cast<int>(std::round(numAmAc * 0.9)); // Correction HH.
			dropletResidues["ATX"] = soluteCharge + static_cast<int>(std::round(numAmAc)) - correction;	// Acetate
			dropletResidues["NXH"] = static_cast<int>(std::round(numAmAc * 0.9));			// Ammonium	
			dropletResidues["NXX"] = static_cast<int>(std::round(numAmAc * 0.1));			// Ammonia
		}
		else {
			std::cerr << "Invalid ESI mode." << std::endl;
			exit(1);
		}

		// Insert new molecules into droplet
		float minSeperation = 0.30f;
		float newRadius = newBoxSize / 2.0f;
		center = { newBoxSize / 2.0f, newBoxSize / 2.0f, newBoxSize / 2.0f };
		insertMolecule("ace.gro",  dropletResidues["ATX"], coordInfo.atoms, coordInfo.coordinates, newRadius, 
			center, minSeperation, proteinCoords, config.droplet_size);
		insertMolecule("aceh.gro", dropletResidues["AHX"], coordInfo.atoms, coordInfo.coordinates, newRadius, 
			center, minSeperation, proteinCoords, config.droplet_size);
		insertMolecule("nh4.gro", dropletResidues["NXH"], coordInfo.atoms, coordInfo.coordinates, newRadius,
			center, minSeperation, proteinCoords, config.droplet_size);
		insertMolecule("nh3.gro",  dropletResidues["NXX"], coordInfo.atoms, coordInfo.coordinates, newRadius, 
			center, minSeperation, proteinCoords, config.droplet_size);

		// Resolvate
		writeGRO("solute.gro", coordInfo.atoms, coordInfo.box_vectors);
		{
			std::ostringstream command;
			std::vector<std::string> inputs = {};
			command << config.gmx_env << " solvate -cp solute.gro -cs tip4p_2005.gro -p temp.top -o solvated.gro";
			auto_gmx_input(command.str(), inputs);
		}
		
		// Expand box for atmosphere
		float finalBoxSize = config.box_size;
		{
			std::ostringstream command;
			command << config.gmx_env << " editconf -f solvated.gro -o expanded.gro -box " <<
				std::to_string(finalBoxSize) << " -c";
			std::vector<std::string> inputs{};
			auto_gmx_input(command.str(), inputs);
		}
		coordInfo = buildCoordInfo("expanded.gro", finalBoxSize);
		
		// Carve droplet
		proteinCoords = getProteinCoords(coordInfo.proteinAtoms);
		std::vector<std::shared_ptr<Atom>> dropletAtoms = carveDroplet(coordInfo, finalBoxSize, config.droplet_size, proteinCoords);
		std::vector<std::array<float, 3>> dropletCoords = extractCoordinates(dropletAtoms);

		// Add in atmosphere (if desired)
		if (config.atm == "yes") {
			coordInfo.atoms = seedAtmosphere(config, dropletAtoms, dropletCoords, coordInfo.box_vectors, newRadius);
		}
		else {
			coordInfo.atoms = dropletAtoms;
		}

		// Correct ordering of atoms to match .top file
		reorderAtoms(coordInfo, topOrder);

		// Write final droplet coord and top files
		createTOP("droplet.top", coordInfo, topOrder);
		writeGRO("droplet.gro", coordInfo.atoms, coordInfo.box_vectors);

		return coordInfo;
	}
}