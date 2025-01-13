#include "Core.h"
#include <vector>

namespace Core {

	// Recenters droplet if drifts more >5nm from COM, prevents interaction with PBC.
	// Modifies coordInfo.atoms directly.
	bool recenterDroplet(CoordInfo& coordInfo, const std::vector<std::array<float, 3>>& proteinCoords) {
		bool recentered = false;

		// Get COM
		std::array<float, 3> com = getCOM(proteinCoords, proteinCoords.size());

		// Get center of box
		std::array<float, 3> boxCenter = { 
			coordInfo.box_vectors[0] / 2.0f, 
			coordInfo.box_vectors[1] / 2.0f,
			coordInfo.box_vectors[2] / 2.0f };


		// If distance is greater than 5nm, recenter droplet
		if (getDistance(com, boxCenter) > 5.0f) {
			recentered = true;

			// Tranlsation vector
			std::array<float, 3> translation = subtractVectors(boxCenter, com);

			// Translate atoms
			for (int idx = 0; idx < coordInfo.atoms.size(); idx++) {
				std::shared_ptr<Atom> atom = coordInfo.atoms[idx];

				// Store which axis atom is outside of box (if any)
				bool outsideLow = false;
				bool outsideHigh = false;
				std::array<bool, 3> outsideLowAxes  = { false, false, false };
				std::array<bool, 3> outsideHighAxes = { false, false, false };

				for (int i = 0; i < 3; i++) {
					atom->coord[i] += translation[i];

					// Determine if gas molecule now outside box
					if (atom->coord[i] < 0.0f) {
						if (atom->res_name == "NNN" || atom->res_name == "OOO") {
							outsideLow = true;
							outsideLowAxes[i] = true;
						}
						else {
							std::cerr << "ERROR: Only gas molecules should be translated" <<
								" across boundaries during recentering." << std::endl;
							std::cerr << atom->toString() << std::endl;
							std::exit(1);
						}
					}
					else if (atom->coord[i] > coordInfo.box_vectors[i]) {
						if (atom->res_name == "NNN" || atom->res_name == "OOO") {
							outsideHigh = true;
							outsideHighAxes[i] = true;
						}
						else {
							std::cerr << "ERROR: Only gas molecules should be translated" <<
								" across boundaries during recentering." << std::endl;
							std::cerr << atom->toString() << std::endl;
							std::exit(1);
						}
					}
				}

				// If gas molecule outside of box, translate entire residue across
				if (outsideLow || outsideHigh) {
					
					// Get both atoms in each gas molecule
					std::vector<std::weak_ptr<Atom>> resAtoms = atom->parent.lock()->atoms;
					int atom_idx;
					if (atom->atom_name == "N1" || atom->atom_name == "O1") {
						atom_idx = idx + 1;
					}
					else if (atom->atom_name == "N2" || atom->atom_name == "O2") {
						atom_idx = idx - 1;
					}
					else {
						std::cerr << "ERROR: Not a gas molecule." << std::endl;
						std::exit(1);
					}

					for (auto& weak_resAtom : resAtoms) {
						auto resAtom = weak_resAtom.lock();
						for (int i = 0; i < 3; i++) {
							if (outsideLowAxes[i]) {
								resAtom->coord[i] += 0.2f;
							}
							else if (outsideHighAxes[i]) {
								resAtom->coord[i] -= 0.2f;
							}
						}
					}

					// Check for steric clashes with reflected gas atoms, if clash, slowly move residue until clash removed
					int iterations = 0;
					while (true) {
						bool clash = false;

						// Find distances with all other gas molecules
						for (std::vector<std::shared_ptr<Atom>>::reverse_iterator riter = coordInfo.atoms.rbegin();
							riter != coordInfo.atoms.rend(); ++riter) {

							// If clash
							auto firstAtom = resAtoms[0].lock();
							if (getDistance(firstAtom->coord, (*riter)->coord) < 1.0f && firstAtom->res_id != (*riter)->res_id) {
								clash = true;

								// Move residue slightly towards center of box incrementally until clash removed
								std::array<float, 3> vector = subtractVectors(firstAtom->coord, boxCenter);
								std::array<float, 3> motion = scaleVector(normalizeVector(vector), 0.1);

								// Move all residue atoms
								for (auto& weak_resAtom : resAtoms) {
									auto resAtom = weak_resAtom.lock();
									for (int i = 0; i < 3; i++) {
										resAtom->coord[i] += motion[i];
									}
								}
							}
						}

						if (!clash) {
							break;
						}
						else if (iterations > 100) {
							std::cerr << "ERROR: Could not remove steric clash between gas molecules during recentering." << std::endl;
							std::exit(1);
						}
						else {
							iterations++;
						}
					}
				}
			}
		}

		return recentered;
	}

}