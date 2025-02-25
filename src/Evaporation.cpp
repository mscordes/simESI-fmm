#include "Core.h"
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <random>

namespace Core {

	// Constants
	float pressure = 1.0f; // 1 atm
	float na = 6.022e23f; // Avogadro's number
	float gasConst = 8.2057366e22f; // Gas constant (nm^3 * atm * K^-1 * mol^-1)

	// Deletes atoms in vector of atoms given vector of indices to delete 
	void deleteAtoms(std::vector<std::shared_ptr<Atom>>& atoms, const std::vector<int>& indicesToDelete) {
		std::vector<int> sortedIndices = indicesToDelete;
		std::sort(sortedIndices.rbegin(), sortedIndices.rend()); // Reverse sort

		// Delete elements at these indices
		for (int index : sortedIndices) {
			if (index >= 0 && index < atoms.size()) { // Ensure the index is valid
				atoms.erase(atoms.begin() + index);
			}
			else {
				std::cerr << "Invalid atom deletion index: " << index << std::endl;
				std::exit(1);
			}
		}
	}

	// Find indicies of molecules that are 'gas phase' ie., not coordinated to protein
	// For example, if atmospheric water present, must distinguis between droplet/atmosperic water
	std::vector<int> isGasPhase(const CoordInfo& coordInfo, const std::vector<std::array<float, 3>>& proteinCarbons, 
		const std::string resType) {

		std::vector<int> vaporIndices;
		float cutoff = 5.0f; // nm - TODO - this may need to be increased for larger droplets
		if (coordInfo.numResidues.at(resType) > 0) {
			for (int i = 0; i < coordInfo.residueMap.at(resType).size(); i++) {
				auto firstAtom = coordInfo.residueMap.at(resType)[i]->atoms[0].lock();
				std::vector<float> dist = norm(proteinCarbons, firstAtom->coord);
				if (*std::min_element(dist.begin(), dist.end()) > cutoff) {
					vaporIndices.push_back(i);
				}
			}
		}
		
		return vaporIndices;
	}

	// Maintains the user inputted gas density of a given solute accounting for evaporation and or solute loss
	static void maintainGasPressure(const Config& config, CoordInfo& coordInfo, const std::vector<std::string>& topOrder,
		const std::string& top_fname, std::vector<std::array<float, 3>>& proteinCarbons, const std::string resType, 
		const int& numGas, const float& fraction, const std::string& resGro, 
		std::unordered_map<std::string, int>& vaporDict, std::vector<int>& evapIndices, 
		std::vector<std::shared_ptr<Atom>>& atomsToInsert) {

		float trueFraction = fraction / 100; // Convert from % to fraction
		int expectedGasNum = static_cast<int>(std::round(trueFraction * numGas)); // Number of this gas to maintain

		// Find solutes that are 'gas phase' ie., not coordinated to protein
		std::vector<int> vaporIndices = isGasPhase(coordInfo, proteinCarbons, resType);

		// Determine if number of solute vapor within cutoff or not, if then delete select number 
		int slop; 
		if (fraction < 0.00001) {
			// If the user defined fraction is 0, delete all evaporated
			slop = 0;
		}
		else {
			slop = static_cast<int>(expectedGasNum * 0.05); // +/- 5% slop
		}
		int delta = expectedGasNum - static_cast<int>(vaporIndices.size());
		if (std::abs(delta) > slop) {

			// Randomly select solutes to delete until expected gas number met
			if (delta < 0) {
				delta = std::abs(delta);
				std::random_device rd;
				std::mt19937 gen(rd());
				std::vector<int> soluteEvapIndices;
				std::sample(vaporIndices.begin(), vaporIndices.end(), std::back_inserter(soluteEvapIndices), delta, gen);

				// Append solutes to be deleted to evapIndices
				for (const auto& solute : soluteEvapIndices) {
					for (const auto& atom : coordInfo.residueMap[resType][solute]->atoms) {
						evapIndices.push_back(atom.lock()->atom_num);
					}
				}

				// Update numResidues and vapor dict
				coordInfo.numResidues[resType] -= delta;
				vaporDict[resType] = static_cast<int>(vaporIndices.size()) - delta;
			}

			// Seed additional solutes if too few
			else {
				int res_num = coordInfo.atoms.back()->res_num;
				for (int i = 0; i < delta; i++) {
					// Randomly select a location for the new molecule
					int attempts = 0; // Insertion attempts
					bool inserted = false;
					while (!inserted) {
						std::random_device rd;
						std::mt19937 gen(rd());
						std::uniform_real_distribution<float> dis(static_cast<float>(0.5f, coordInfo.box_vectors[0] - 0.5f));

						// Generate random coordinates
						std::array<float, 3> newLoc = { dis(gen), dis(gen), dis(gen) };

						// New molecule atoms
						std::vector<std::shared_ptr<Atom>> molecAtoms = readGROatoms(resGro);

						// Molecule coordinates
						for (auto& atom : molecAtoms) {
							for (int j = 0; j < 3; j++) {
								atom->coord[j] += newLoc[j];
							}
						}

						// Check if the molecule is too close to another molecule
						bool tooClose = false;
						float minSepSquared = 0.5f * 0.5f;
						for (const auto& atom : molecAtoms) {
							std::array<float, 3> newCoord = atom->coord;
							for (const auto& coord : coordInfo.coordinates) {
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
						if (!tooClose) {
							res_num++;
							for (const auto& atom : molecAtoms) {
								coordInfo.coordinates.push_back(atom->coord);
								atom->res_num = res_num;
								atom->res_id = atom->res_name + std::to_string(res_num);
								atom->velocity = sampleMaxwell(1, config.gas_temp, atom->element)[0];
								atomsToInsert.push_back(atom);
							}
							inserted = true;
						}

						// If too many attempts, break
						attempts++;
						if (attempts > 10000) {
							std::cerr << "ERROR: Could not insert all molecules during solute vapor seeding." << std::endl;
							std::exit(1);
						}
					}
				}

				// Update numResidues and vapor dict
				coordInfo.numResidues[resType] += delta;
				vaporDict[resType] = static_cast<int>(vaporIndices.size()) + delta;
			}
		}

		// Else if delta less than slop, num vapor is vapor indices
		else {
			vaporDict[resType] = static_cast<int>(vaporIndices.size());
		}
	}

	// Find and remove evaporated waters and solutes while maintaing constant (user defined) humidity level
	bool removeEvaporated(const Config& config, CoordInfo& coordInfo, const std::vector<std::string>& topOrder, 
		const std::string& top_fname, std::vector<std::array<float, 3>>& proteinCarbons, std::unordered_map<std::string, int>& vaporDict) {

		bool evap = false; // True if any evaporated molecules deleted

		// Get coordinates of protein carbons
		for (const auto& [mononmer, monAtoms] : coordInfo.proteinAtoms) {
			for (const auto& atom : monAtoms) {
				if (atom->atom_name == "C") {
					proteinCarbons.push_back(atom->coord);
				}
			}
		}

		// Find number of total gas present 
		float boxVolume = coordInfo.box_vectors[0] * coordInfo.box_vectors[1] * coordInfo.box_vectors[2];
		float numDensity = (pressure * na) / (gasConst * config.gas_temp);
		int numGas = static_cast<int>(std::round(numDensity * boxVolume));
		
		// Fraction of each gas type
		std::unordered_map<std::string, float> fractionDict = {
			{ "SOL", config.water_vapor }, { "ATX", config.ace_vapor }, { "AHX", config.ach_vapor }, 
			{ "NXH", config.nh4_vapor }, { "NXX", config.nh3_vapor }, { "OHX", 0.0f }, { "HHO", 0.0f }
		};

		// Number of N2 and O2
		float summed_fraction = 0.0f;
		for (const auto& pair : fractionDict) {
			summed_fraction += pair.second;
		}
		fractionDict["NNN"] = (100.0f - summed_fraction) * 0.786f; // 78.6% N2 in atmosphere
		fractionDict["OOO"] = (100.0f - summed_fraction) * 0.214f; // 21.4% O2 in atmosphere

		// Names of gro files to extract coordinate file information from (if need to seed additional)
		std::unordered_map<std::string, std::string> groDict = {
			{ "SOL", "sol.gro"}, { "ATX", "ace.gro"}, { "AHX", "aceh.gro"},
			{ "NXH", "nh4.gro"}, { "NXX", "nh3.gro"}, { "OHX", "oh.gro"}, { "HHO", "h3o.gro" }, 
			{ "NNN", "n2.gro"}, { "OOO", "o2.gro"}
		};

		// Dictionary to keep track of number of each gas type in 'gas' phase
		vaporDict = {
			{ "SOL", 0 }, { "ATX", 0 }, { "AHX", 0 }, { "NXH", 0 }, { "NXX", 0 }, { "OHX", 0 }, { "HHO", 0 }, 
			{ "NNN", 0 }, { "OOO", 0 }
		};

		// Actual insertion of deletion of gas molecules to maintain user inputted pressure of each gas
		std::vector<std::shared_ptr<Atom>> atomsToInsert;
		std::vector<int> evapIndices; 
		for (const auto& [ resType, fraction ] : fractionDict) {
			maintainGasPressure(config, coordInfo, topOrder, top_fname, proteinCarbons, resType, numGas, fraction,
				groDict[resType], vaporDict, evapIndices, atomsToInsert);
		}

		// Remove evaporated atoms and insert new atoms (if necessary)
		if (!evapIndices.empty() || !atomsToInsert.empty()) {
			evap = true;

			// Delete first
			if (!evapIndices.empty()) {
				deleteAtoms(coordInfo.atoms, evapIndices);

				// Ensure correct atom numbering post evaporation
				for (int i = 0; i < coordInfo.atoms.size(); i++) {
					coordInfo.atoms[i]->atom_num = i;
				}
			}

			// Insert new atoms (ensuring correct ordering)
			if (!atomsToInsert.empty()) {
				coordInfo.atoms.insert(coordInfo.atoms.end(), atomsToInsert.begin(), atomsToInsert.end());
				reorderAtoms(coordInfo, topOrder);
			}

			// Final upate of .top file
			createTOP(top_fname, coordInfo, topOrder);
		}

		return evap;
	}

}