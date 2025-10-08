#include "Recentering.h"
#include "Atom.h"
#include "Atmosphere.h"
#include "CoordinateManip.h"
#include "MathUtils.h"

void recenter(
	State& state,
	unordered_set<shared_ptr<Residue>>& gWaters,
	RunFlags& flags,
	double cutoff) 
{
	static const double boxD = Parameters::Get().getBoxSize();
	static double half = Parameters::Get().getBoxSize() / 2;
	static const Vec3D bCenter = { half, half, half };

	const auto pCoords = getProteinCoords(state);
	const Vec3D pCenter = getCenter(pCoords);

	if (pCenter.distance(bCenter) < cutoff) return;
	flags.createRunFile = true;

	// Begin recentering
	Vec3D translation = bCenter - pCenter;
	for (const auto& [rType, residues] : state.residueSet)
		for (const auto& r : residues)
			for (const auto& a : r->atoms)
				a->coord = a->coord + translation;

	// Check if gas molecules are outside of box, if so, reflect across boundary
	static const auto gasTypes = getGasNums();
	for (const auto& [gType, _] : gasTypes) {
		auto it = state.residueSet.find(gType);
		if (it == state.residueSet.end()) continue;

		unordered_set<shared_ptr<Residue>> gases;
		if (gType != "SOL") gases = it->second;
		else				gases = gWaters;
		if (gases.empty()) continue;

		for (const auto& g : gases) {
			const auto& atoms = g->atoms;

			Vec3D gTranslation = { 0.0, 0.0, 0.0 };
			bool cx = false, cy = false, cz = false;
			for (const auto& a : atoms) {
				auto& c = a->coord;
				if (!cx && c.x < 0) {
					cx = true;
					gTranslation.x = gTranslation.x + boxD;
				}
				if (!cy && c.y < 0) {
					cy = true;
					gTranslation.y = gTranslation.y + boxD;
				}
				if (!cz && c.z < 0) {
					cz = true;
					gTranslation.z = gTranslation.z + boxD;
				}

				if (!cx && c.x > boxD) {
					cx = true;
					gTranslation.x = gTranslation.x - boxD;
				}
				if (!cy && c.y > boxD) {
					cy = true;
					gTranslation.y = gTranslation.y - boxD;
				}
				if (!cz && c.z > boxD) {
					cz = true;
					gTranslation.z = gTranslation.z - boxD;
				}
			}

			if (cx || cy || cz) {
				for (const auto& a : atoms)
					a->coord = a->coord + gTranslation;

				// Now check if gas molecules are 'split' across boundaries, if so, move towards center of box
				int iterations = 0;
				while (true) {
					++iterations;
					if (iterations > 1000)
						throw runtime_error("recenter: Could not move gas molecule to free position during recentering.");

					bool outside = false;
					bool clash = false;

					for (const auto& a : atoms) {
						auto& c = a->coord;

						// Whether any atoms still outside box
						for (const auto& i : { c.x, c.y, c.z }) {
							if (i < 0 || i > boxD) {
								outside = true;
								goto endloop;
							}
						}

						// Check for clashes with other gas molecules
						for (const auto& ngType : { "NNN", "OOO" }) {
							auto nIt = state.residueSet.find(ngType);
							if (nIt == state.residueSet.end()) continue;
							const auto& ngases = nIt->second;

							for (const auto& ng : ngases) {
								if (ng.get() == g.get()) continue; // skip self
								for (const auto& na : ng->atoms) {
									if (c.distance_sq(na->coord) < 1.0) {
										clash = true;
										goto endloop;
									}
								}
							}
						}
					}
					endloop:

					// Move residue slightly towards center of box incrementally until clash removed
					Vec3D nTranslation = (bCenter - atoms[0]->coord).normalize() * 0.1;
					for (const auto& a : atoms)
						a->coord = a->coord + nTranslation;

					if (!outside && !clash) break;
				}
			}
		}
	}
}