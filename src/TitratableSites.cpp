#include "Core.h"
#include <algorithm>
#include <memory>

namespace Core {

    // pdb2gmx takes residues types in specific order, exclude termini for now
    const std::vector<std::string> proteinResNames = { "LYS", "ARG", "ASP", "GLU", "HIS" };

    // Titrable hydrogens of each titrable amino acid (CHARMM36)
    const std::unordered_map<std::string, std::vector<std::string>> proteinDonorsMap = {
        { "LYS", {"HZ1", "HZ2", "HZ3"} }, { "ARG", {"HH11", "HH12", "HH21", "HH22"} },
        { "ASP", {"HD2"} }, { "GLU", {"HE2"} }, { "HIS", {"HD1", "HE2"} },
        {"NTERM", {"H1", "H2", "H3"}},  {"CTERM", {"HT2"}}
    };

    // Titratatable hydrogen acceptor names of each titrable amino acid (CHARMM36)
    const std::unordered_map<std::string, std::vector<std::string>> proteinAcceptorsMap = {
        { "LYS", {"NZ"} }, { "ARG", {"NH1"} }, { "ASP", {"OD1", "OD2"}}, {"GLU", {"OE1", "OE2"}}, {"HIS", {"NE2", "ND1"}},
        {"NTERM", {"N"}},  {"CTERM", {"OT1", "OT2"}}
    };

    // Titratatable hydrogen names for non-protein residues
    const std::unordered_map<std::string, std::vector<std::string>> donorsMap = {
        {"HHO", {"HW1", "HW2", "HW3"}}, {"AHX", {"HO1"}}, {"NXH", {"HZ1", "HZ2", "HZ3", "HZ4"}}
    };

    // Titratatable hydrogen acceptor names for non-protein residues
    const std::unordered_map<std::string, std::vector<std::string>> acceptorsMap = {
        { "OHX", {"O1"}}, {"ATX", {"O1", "O2"}}, {"NXX", {"N1"}}
    };

    // Non-protein pKa values
	const std::unordered_map<std::string, float> nonProteinPkas = {
        {"HHO", 0.0f}, {"OHX", 14.0f}, {"ATX", 4.76f}, {"AHX", 4.76f}, {"NXX", 9.25f}, {"NXH", 9.25f}
	};

    // Gas phase basicities (kJ/mol)
    const std::unordered_map<std::string, float> gasPhaseBasicities = {
        {"HIS", 935.54f}, {"LYS", 884.08f}, {"ARG", 983.24f}, {"GLU", 1424.23f}, {"ASP", 1428.42f}, {"NTERM", 850.90f},
        {"CTERM", 1400.38f}, {"OHX", 1605.40f}, {"HHO", 659.82f}, {"ATX", 1428.42f}, {"AHX", 1428.42f}, {"NXX", 818.81f}, {"NXH", 818.81f}
    };

	// Determine protonation state of a protein resiude type and store either hydrogen donors or acceptors
    static void proteinTitAtoms(const std::string& resType, const std::vector<std::shared_ptr<Residue>>& residues,
        std::vector<std::shared_ptr<Atom>>& donors, std::vector<std::shared_ptr<Atom>>& acceptors) {

        for (const auto& residue : residues) {
            std::vector<std::shared_ptr<Atom>> potentialDonors;
            std::vector<std::shared_ptr<Atom>> potentialAcceptors;
            std::vector<std::string> hydrogenNames = proteinDonorsMap.at(resType);
            std::vector<std::string > acceptorNames = proteinAcceptorsMap.at(resType);

            for (const auto& weak_atom : residue->atoms) {
                auto atom = weak_atom.lock();
                if (std::find(hydrogenNames.begin(), hydrogenNames.end(), atom->atom_name) != hydrogenNames.end()) {
                    potentialDonors.push_back(atom);
                }
                else if (std::find(acceptorNames.begin(), acceptorNames.end(), atom->atom_name) != acceptorNames.end()) {
                    potentialAcceptors.push_back(atom);
                }
            }

            // If all hydrogens are present, add all to donors
            if (potentialDonors.size() == hydrogenNames.size()) {
                for (const auto& donor : potentialDonors) {
                    donors.push_back(donor);
                }
            }
            else {
                // HIS can have 3 protonation states, so parse here
                if (residue->res_name == "HIS") {
                    for (const auto& weak_atom : residue->atoms) {
                        auto atom = weak_atom.lock();

                        // HISD
                        if (atom->atom_name == "HD1") {
                            for (const auto& acceptor : potentialAcceptors) {
                                if (acceptor->atom_name == "NE2") {
                                    acceptors.push_back(acceptor);
                                }
                            }
                            break;
                        }
                        // HISE
                        else if (atom->atom_name == "HE2") {
                            for (const auto& acceptor : potentialAcceptors) {
                                if (acceptor->atom_name == "ND1") {
                                    acceptors.push_back(acceptor);
                                }
                            }
                            break;
                        }
                    }
                }
                else {
                    for (const auto& acceptor : potentialAcceptors) {
                        acceptors.push_back(acceptor);
                    }
                }
            }
        }
    }

    // Finds tit atoms of non-protein residues (prot state determined by inputted map and atoms list)
    static void nonProteinTitAtoms(std::vector<std::shared_ptr<Atom>>& atoms, const CoordInfo& coordInfo, 
        const std::unordered_map<std::string, std::vector<std::string>>& map) {

        for (auto const& [resType, atomNames] : map) {
            for (const auto& residue : coordInfo.residueMap.at(resType)) {
                for (const auto& weak_atom : residue->atoms) {
                    auto atom = weak_atom.lock();
                    if (std::find(atomNames.begin(), atomNames.end(), atom->atom_name) != atomNames.end()) {
                        atoms.push_back(atom);
                    }
                }
            }
        }
    }

	// Pin gas phase basicities to atoms
    static std::vector<float> pinGPB(const std::vector<std::shared_ptr<Atom>>& atoms) {
        std::vector<float> gpbs;

        for (const auto& atom : atoms) {
            if (gasPhaseBasicities.find(atom->res_name) != gasPhaseBasicities.end()) {
                gpbs.push_back(gasPhaseBasicities.at(atom->res_name));
			}
			else {
                // Termini can be present in any amino acid, so check here
				if (std::find(proteinDonorsMap.at("NTERM").begin(), proteinDonorsMap.at("NTERM").end(), atom->atom_name) != proteinDonorsMap.at("NTERM").end()
                    || std::find(proteinAcceptorsMap.at("NTERM").begin(), proteinAcceptorsMap.at("NTERM").end(), atom->atom_name) != proteinAcceptorsMap.at("NTERM").end()) {
					gpbs.push_back(gasPhaseBasicities.at("NTERM"));
				}
                else if (std::find(proteinDonorsMap.at("CTERM").begin(), proteinDonorsMap.at("CTERM").end(), atom->atom_name) != proteinDonorsMap.at("CTERM").end()
                    || std::find(proteinAcceptorsMap.at("CTERM").begin(), proteinAcceptorsMap.at("CTERM").end(), atom->atom_name) != proteinAcceptorsMap.at("CTERM").end()) {
                    gpbs.push_back(gasPhaseBasicities.at("CTERM"));
                }
                else {
					std::cerr << "Error: Could not find gas phase basicity for atom " << atom->toString() << std::endl;
                    exit(1);
                }
			}
        }

        return gpbs;
    }

    // Pin pKa's to atoms
    static std::vector<float> pinPKA(const std::vector<std::shared_ptr<Atom>>& atoms, const std::unordered_map<int, float>& pkaMap) {
        std::vector<float> pkas;

        for (const auto& atom : atoms) {

			// Non-protein pKa's
			if (nonProteinPkas.find(atom->res_name) != nonProteinPkas.end()) {
				pkas.push_back(nonProteinPkas.at(atom->res_name));
			}

            // For protein, pin individual pKa's as computed via PROPKA3
			else {
				if (pkaMap.find(atom->res_num) == pkaMap.end()) {
					std::cerr << "Error: Could not find pKa for atom " << atom->toString() << std::endl;
					exit(1);
				}
				else {
					pkas.push_back(pkaMap.at(atom->res_num));
				}
			}
        }

        return pkas;
    }

	// Bin all (non-water) proton donors/acceptors into TitratableAtom objects 
    TitratableSites getTitratableSites(const CoordInfo& coordInfo, const std::vector<int>&clusterIDs,
        const std::unordered_map<int, Cluster>& clusters, const std::unordered_map<int, float>&pkaMap) {

        // Get titratable atoms first
        std::vector<std::shared_ptr<Atom>> donorAtoms;
        std::vector<std::shared_ptr<Atom>> acceptorAtoms;

        // Start with amino acids. Determine protonation state, and add bin hydrogen donors/acceptors
        for (const auto& resType : proteinResNames) {
			std::vector<std::shared_ptr<Residue>> residues = coordInfo.residueMap.at(resType);
			proteinTitAtoms(resType, residues, donorAtoms, acceptorAtoms);
		}

		// Don't forget to add N and C termini
        for (const auto& monomer : coordInfo.proteins) {
            std::vector<std::shared_ptr<Residue>> nterm = { monomer->residues[0].lock()};
            proteinTitAtoms("NTERM", nterm, donorAtoms, acceptorAtoms);

            std::vector<std::shared_ptr<Residue>> cterm = { monomer->residues.back().lock()};
            proteinTitAtoms("CTERM", cterm, donorAtoms, acceptorAtoms);
        }

		// Now add non-protein, non-water titratable hydrogens
		nonProteinTitAtoms(donorAtoms, coordInfo, donorsMap);
		nonProteinTitAtoms(acceptorAtoms, coordInfo, acceptorsMap);

		// Get coordinates of donor/acceptor atoms
		std::vector<std::array<float, 3>> donorCoords = extractCoordinates(donorAtoms);
		std::vector<std::array<float, 3>> acceptorCoords = extractCoordinates(acceptorAtoms);

        /* Now that we have atoms, determine ID of cluster solvating atom, as well as coordinated waters
        (ie, the number of waters in the solvating cluster) */
		std::unordered_map<int, int> clusterWaters = findClusterWaters(clusterIDs);
		std::vector<int> donorClusterIDs = pinCluster(donorAtoms, coordInfo.waterOCoords, clusterIDs, clusterWaters, true);
		std::vector<int> donorNearWaters = pinClusterWaters(donorClusterIDs, clusters);
		std::vector<int> acceptorClusterIDs = pinCluster(acceptorAtoms, coordInfo.waterHCoords, clusterIDs, clusterWaters, false);
        std::vector<int> acceptorNearWaters = pinClusterWaters(acceptorClusterIDs, clusters);

        // Determine pKa and gas phase basicities (GPB's)
		std::vector<float> donorPkas = pinPKA(donorAtoms, pkaMap);
		std::vector<float> donorGpbs = pinGPB(donorAtoms);
		std::vector<float> acceptorPkas = pinPKA(acceptorAtoms, pkaMap);
		std::vector<float> acceptorGpbs = pinGPB(acceptorAtoms);

		// Bin into TitratableAtom objects, and collimate into TitratableSites object
		std::vector<TitratableAtom> donors = {};
		for (size_t i = 0; i < donorAtoms.size(); i++) {
			donors.push_back(TitratableAtom(donorAtoms[i], donorCoords[i], donorNearWaters[i], donorClusterIDs[i], donorPkas[i], donorGpbs[i]));
		}
        std::vector<TitratableAtom> acceptors = {};
        for (size_t i = 0; i < acceptorAtoms.size(); i++) {
            acceptors.push_back(TitratableAtom(acceptorAtoms[i], acceptorCoords[i], acceptorNearWaters[i], acceptorClusterIDs[i], acceptorPkas[i], acceptorGpbs[i]));
        }

		return TitratableSites(donors, acceptors);
	}
}