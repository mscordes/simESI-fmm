#include "Core.h"
#include <algorithm>
#include <unordered_map>

namespace Core {

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

	// Find and remove evaporated waters and solutes
	bool removeEvaporated(CoordInfo& coordInfo, const std::vector<std::string>& topOrder, const std::string& top_fname, 
		std::vector<std::array<float, 3>>& proteinCarbons) {
		bool evap = false; // True if any evaporated molecules deleted
		float cutoff = 10.0; //nm

		// Get coordinates of protein carbons
		for (const auto& [mononmer, monAtoms] : coordInfo.proteinAtoms) {
			for (const auto& atom : monAtoms) {
				if (atom->atom_name == "C") {
					proteinCarbons.push_back(atom->coord);
				}
			}
		}

		std::vector<int> evapIndices;
		for (const auto& resType : topOrder) {
			if (resType != "NNN" && resType != "OOO" && coordInfo.numResidues[resType] > 0) {
				const auto& residues = coordInfo.residueMap[resType];
				for (const auto& residue : residues) {
					auto firstAtom = residue->atoms[0].lock();
					std::vector<float> dist = norm(proteinCarbons, firstAtom->coord);
					if (*std::min_element(dist.begin(), dist.end()) > cutoff) {

						// Update numResidues
						coordInfo.numResidues[resType] -= 1;

						for (const auto& atom : residue->atoms) {
							evapIndices.push_back(atom.lock()->atom_num);
						}
					}
				}
			}
		}

		// Remove evaporated waters and solutes
		if (!evapIndices.empty()) {
			evap = true;

			// Update atoms and .top, dont need to rebuild coordInfo since we are not using it after this
			deleteAtoms(coordInfo.atoms, evapIndices);
			createTOP(top_fname, coordInfo, topOrder);
		}

		return evap;
	}
}