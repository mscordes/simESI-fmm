#include "Core.h"
#include <vector>

namespace Core {

	// Determines if a residue is outside the boundaries of the simulation box, if so, by which axes
	static void isOutside(const CoordInfo& coordInfo, const std::shared_ptr<Residue>& residue, bool& outsideLow, bool& outsideHigh, 
		std::array<bool, 3>& outsideLowAxes, std::array<bool, 3>& outsideHighAxes) {

		std::vector<std::weak_ptr<Atom>> resAtoms = residue->atoms;
		for (auto& weakAtom : resAtoms) {
			auto atom = weakAtom.lock();

			for (int i = 0; i < 3; i++) {
				if (atom->coord[i] < 0.0f) {
					if (atom->res_name == "NNN" || atom->res_name == "OOO" || atom->res_name == "SOL") {
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
					if (atom->res_name == "NNN" || atom->res_name == "OOO" || atom->res_name == "SOL") {
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
		}
	}

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
				for (int i = 0; i < 3; i++) {
					atom->coord[i] += translation[i];
				}
			}

			// Check if gas molecules are outside of box, if so, reflect across boundary
			for (auto& residue : coordInfo.residues) {
				bool outsideLow = false;
				bool outsideHigh = false;
				std::array<bool, 3> outsideLowAxes = { false, false, false };
				std::array<bool, 3> outsideHighAxes = { false, false, false };
				isOutside(coordInfo, residue, outsideLow, outsideHigh, outsideLowAxes, outsideHighAxes);

				// Reflect based on which axes are outside
				if (outsideLow || outsideHigh) {

					// Reflect atoms across boundary of box
					std::vector<std::weak_ptr<Atom>> resAtoms = residue->atoms;
					for (auto& weak_resAtom : resAtoms) {
						auto resAtom = weak_resAtom.lock();
						for (int i = 0; i < 3; i++) {
							if (outsideLowAxes[i]) {
								resAtom->coord[i] += coordInfo.box_vectors[i];
							}
							else if (outsideHighAxes[i]) {
								resAtom->coord[i] -= coordInfo.box_vectors[i];
							}
						}
					}
				}
			}

			// Now check if gas molecules are 'split' across boundaries, if so, move towards center of box
			for (auto& residue : coordInfo.residues) {
				bool outsideLow = false;
				bool outsideHigh = false;
				std::array<bool, 3> outsideLowAxes = { false, false, false };
				std::array<bool, 3> outsideHighAxes = { false, false, false };
				isOutside(coordInfo, residue, outsideLow, outsideHigh, outsideLowAxes, outsideHighAxes);

				// Determine which axes are split outside
				if (outsideLow || outsideHigh) {

					// Give 'push' until atoms are fully in box (0.5nm increments)
					std::vector<std::weak_ptr<Atom>> resAtoms = residue->atoms;
					for (auto& weak_resAtom : resAtoms) {
						auto resAtom = weak_resAtom.lock();
						for (int i = 0; i < 3; i++) {
							if (outsideLowAxes[i]) {
								resAtom->coord[i] += 0.5f;
							}
							else if (outsideHighAxes[i]) {
								resAtom->coord[i] -= 0.5f;
							}
						}
					}
									
					// Check if pushed gas atoms have clash, if so, slowly move residue towards center until clash removed
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