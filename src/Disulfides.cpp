#include "Core.h"
#include <cmath>

namespace Core {

	/*Annoying bug where will pdb2gmx will not recognize disulfide bond if bond extends beyond 2.0 ± 0.2 Å.
    Not only do we want to keep disulfides intact, but simESI requires every atom be accounted for so this 
    bug actually breaks the simulation as well.*/
    std::vector<std::shared_ptr<Atom>> fixDisulfides(const std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>>& residueMap,
        std::vector<std::shared_ptr<Atom>>& atoms) {

		if (residueMap.find("CYS") == residueMap.end()) {
			return atoms;
		}
        else {
            // Get indices of sulfur atoms
            std::vector<std::shared_ptr<Atom>> ss_indices {};
            for (const auto& cys : residueMap.at("CYS")) {
				for (const auto& atom : cys->atoms) {
					if (atom.lock()->atom_name == "SG") {
						ss_indices.push_back(atom.lock());
					}
				}
            }

			// Find sulfur atoms that are within 2.5 Å
			std::vector<std::pair<Atom, Atom>> ss_pairs {};
            std::vector<Atom> void_atoms {};
			for (const auto& atom1 : ss_indices) {
				for (const auto& atom2 : ss_indices) {
					if (atom1->res_id == atom2->res_id) continue; // Avoid duplicate pairs and self-pairing
					float dist = getDistance(atom1->coord, atom2->coord);
					if (dist < 0.25f) {
						if (dist < 0.19f || dist > 0.21f) {
							std::array<float, 3> com = getCOM({ atom1->coord, atom2->coord }, 2);
							atom1->coord = setBondLength(com, atom1->coord, 0.1025f);
							atom2->coord = setBondLength(com, atom2->coord, 0.1025f);
						}
					}
				}
			}

			return atoms;
        }
    }
}