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

	// Get number of each gas type to seed
	static int getNumGas(std::vector<int>& composition, const int& totGas, const float& fraction, const int& mapVal) {\
		int num = static_cast<int>((fraction / 100) * totGas);
		for (int i = 0; i < num; i++) {
			composition.push_back(mapVal);
		}
		return num;
	}

	// Seed in atmosphere
	std::vector<std::shared_ptr<Atom>> seedAtmosphere(const Config& config, std::vector<std::shared_ptr<Atom>>& atoms, 
		std::vector<std::array<float, 3>>& coords, const std::array<float, 3>& boxVectors, const float& dropletRadius) {

		// Find number of gas to seed
		float pressure = 1.0f; // 1 atm
		float na = 6.022e23f; // Avogadro's number
		float gasConst = 8.2057366e22f; // Gas constant (nm^3 * atm * K^-1 * mol^-1)
		float boxVolume = boxVectors[0] * boxVectors[1] * boxVectors[2];
		float numDensity = (pressure * na) / (gasConst * config.gas_temp);
		int numGas = static_cast<int>(std::round(numDensity * boxVolume));

		/* For initial gas positions, seed in grid then randomly pertube allowing only some displacement
		because need ints, rounding will lead to less than gas_num being seeded so scale up until threshold met */
		float scale = 1.0f;
		int numPoints = static_cast<int>(std::round(std::pow(numGas, 1.0f / 3.0f)));
		while (true) {
			if (numPoints * numPoints * numPoints < numGas) {
				scale += 0.1f;
				numPoints = static_cast<int>(std::round(std::pow(numGas / scale, 1.0f / 3.0f)));
			}
			else {
				break;
			}
		}

		// Create grid
		float startPoint = 0.3f; // Buffer from edge
		float endPoint = boxVectors[0] - 0.3f;
		float step = (endPoint - startPoint) / numPoints;

		std::vector<std::array<float, 3>> gasCoords;
		for (float x = startPoint; x < endPoint; x += step) {
			for (float y = startPoint; y < endPoint; y += step) {
				for (float z = startPoint; z < endPoint; z += step) {
					std::array<float, 3> newCoord = { x, y, z };
					gasCoords.push_back(newCoord);
				}
			}
		}

		// Amount gases can move during pertubation
		float minSpace = 0.6f; // (nm), closest distance gases can be, ~2 N2 widths
		float maxMove = (step - minSpace) / 2.0f;

		// Perturb gas positions 
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> dis(-maxMove, maxMove);
		for (auto& coord : gasCoords) {
			std::array<float, 3> perturb = { dis(gen), dis(gen), dis(gen) };
			std::array<float, 3> newCoord = { coord[0] + perturb[0], coord[1] + perturb[1], coord[2] + perturb[2] };
			coord = newCoord;
		}

		// Remove coordinates that overlap with droplet
		std::vector<int> removeIndices;
		float cutoff = dropletRadius * 1.5f;
		for (int i = 0; i < gasCoords.size(); i++) {
			if (getDistance(gasCoords[i], { boxVectors[0] / 2.0f, boxVectors[1] / 2.0f, boxVectors[2] / 2.0f }) < cutoff) {
				removeIndices.push_back(i);
			}
		}

		// Because of scaling, need to remove some points
		int toRemove = static_cast<int>(gasCoords.size()) - numGas - static_cast<int>(removeIndices.size());
		if (toRemove > 0) {
			for (int i = 0; i < toRemove; i++) {
				while (true) {
					std::uniform_int_distribution<int> dis(0, static_cast<int>(gasCoords.size()) - 1);
					int index = dis(gen);
					if (std::find(removeIndices.begin(), removeIndices.end(), index) == removeIndices.end()) {
						removeIndices.push_back(index);
						break;
					}
				}
			}
		}

		// Remove duplicates in removeIndices
		std::set<int> s(removeIndices.begin(), removeIndices.end());
		removeIndices.assign(s.begin(), s.end());

		// Remove indices in reverse order
		std::sort(removeIndices.begin(), removeIndices.end(), std::greater<>());
		for (const auto& index : removeIndices) {
			gasCoords[index] = gasCoords.back();
			gasCoords.pop_back();
		}
		numGas = gasCoords.size();

		// Determine how many of each gas type to seed including water
		std::vector<int> composition;
		int num_h2o = getNumGas(composition, numGas, config.water_vapor, 2); 
		int num_ace = getNumGas(composition, numGas, config.ace_vapor, 3); 
		int num_ach = getNumGas(composition, numGas, config.ach_vapor, 4); 
		int num_nh4 = getNumGas(composition, numGas, config.nh4_vapor, 5); 
		int num_nh3 = getNumGas(composition, numGas, config.nh3_vapor, 6); 
		int num_n2_o2 = numGas - num_h2o - num_ace - num_ach - num_nh4 - num_nh3;
		int num_n2 = static_cast<int>(num_n2_o2 * 0.786f); // 78.6% N2 in atmosphere
		for (int i = 0; i < num_n2; i++) {
			composition.push_back(0);
		}
		int num_o2 = num_n2_o2 - num_n2; // 21.4% O2 in atmosphere
		for (int i = 0; i < num_o2; i++) {
			composition.push_back(1);
		}

		// Final shuffle
		std::shuffle(composition.begin(), composition.end(), gen);

		// .gro names of each molecule type
		std::unordered_map<int, std::string> groDict = { 
			{ 2, "sol.gro" },
			{ 3, "ace.gro" }, 
			{ 4, "aceh.gro" },
			{ 5, "nh4.gro" },
			{ 6, "nh3.gro" },
		};

		// Seed gas
		std::string res_id;
		int res_num = atoms.back()->res_num;
		std::string res_name;
		std::string atom_name;
		std::string element;
		for (int i = 0; i < numGas; ++i) {
			res_num++;

			if (composition[i] == 0) {
				res_name = "NNN";
				res_id = res_name + std::to_string(res_num);
				element = "N";
				for (int j = 0; j < 2; j++) {
					atom_name = "N" + std::to_string(j + 1);

					if (j == 1) {
						gasCoords[i][0] += 0.10977f; //N2 bond length
					}

					atoms.push_back(std::make_shared<Atom>
						(nullptr, res_id, res_num, res_name, atom_name, 0, gasCoords[i], std::array<float, 3>{ 0.0f, 0.0f, 0.0f }, element, "Z"));
				}
			}
			else if (composition[i] == 1) {
				res_name = "OOO";
				res_id = res_name + std::to_string(res_num);
				element = "O";
				for (int j = 0; j < 2; j++) {
					atom_name = "O" + std::to_string(j + 1);

					if (j == 1) {
						gasCoords[i][0] += 0.12074f; //O2 bond length
					}

					atoms.push_back(std::make_shared<Atom>
						(nullptr, res_id, res_num, res_name, atom_name, 0, gasCoords[i], std::array<float, 3>{ 0.0f, 0.0f, 0.0f }, element, "Z"));
				}
			}
			else {
				// Seed other non N2/O2 gases
				std::string groFile = groDict[composition[i]];
				std::vector<std::shared_ptr<Atom>> molecAtoms = readGROatoms(groFile);
				std::vector<std::array<float, 3>> molecCoords = extractCoordinates(molecAtoms);
				for (std::shared_ptr<Atom> atom : molecAtoms) {

					for (int k = 0; k < 3; k++) {
						atom->coord[k] += gasCoords[i][k];
					}
					atom->res_num = res_num;
					atom->res_id = atom->res_name + std::to_string(res_num);
					atoms.push_back(atom);
				}
			}
		}

		return atoms;
	}
}