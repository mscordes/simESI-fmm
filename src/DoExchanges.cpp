#include "Core.h"
#include <algorithm>
#include <unordered_map>

/* Facilitates exchanges accepted by foundPairs(). */

namespace Core {

	// Swaps iter of two atoms in a vector
	static void swapAtoms(std::vector<std::shared_ptr<Atom>>& atoms, int i, int j) {
		if (i >= atoms.size() || j >= atoms.size()) {
			std::cerr << "Invalid indices for atom swapping: " << i << ", " << j << std::endl;
			exit(1);
		}
		std::swap(atoms[i], atoms[j]);
	}

	// Swaps coordinates and velocites of two atoms. Modifies the objects in place.
	static void swapCoordVel(std::shared_ptr<Atom>& atom1, std::shared_ptr<Atom>& atom2) {
		std::array<float, 3> tempCoord = atom1->coord;
		std::array<float, 3> tempVel = atom1->velocity;
		atom1->coord = atom2->coord;
		atom1->velocity = atom2->velocity;
		atom2->coord = tempCoord;
		atom2->velocity = tempVel;
	}

	// Inserts new protein protons while preserving iter
	static void insertAtoms(std::vector<std::shared_ptr<Atom>>& atoms, const std::vector<std::shared_ptr<Atom>>& newHydrogens,
		const std::unordered_map<std::string, std::string>& proteinAcceptingMap) {

		// Residue and atoms accepting protons, map keys are accepting residue name + accepting atom name
		std::vector<std::string> parentResidues;
		std::vector<std::string> parentAtoms;
		std::unordered_map<std::string, std::shared_ptr<Atom>> parentAtomsMap;
		for (const auto& hydrogen : newHydrogens) {
			// Accpeting residue
			parentResidues.push_back(hydrogen->res_id);

			// Begin finding accepting atom name for each new hydrogen
			std::string parentAtomName;

			// Check for termini
			if ((hydrogen->atom_name == "H3" && proteinAcceptingMap.find(hydrogen->res_id) != proteinAcceptingMap.end())
				|| (hydrogen->atom_name == "HT2" && proteinAcceptingMap.find(hydrogen->res_id) != proteinAcceptingMap.end())) {
				parentAtomName = proteinAcceptingMap.at(hydrogen->res_id);
			}
			// HIS, take care for 2 protonation states
			else if (hydrogen->res_name.find("HIS") != std::string::npos) {
				// HISE
				if (hydrogen->atom_name == "HE2") {
					parentAtomName = "NE2";
				}
				// HISD
				else if (hydrogen->atom_name == "HD1") {
					parentAtomName = "ND1";
				}
				else {
					std::cerr << "ERROR: HIS donor H name not recognized." << std::endl;
					std::cerr << hydrogen->toString() << std::endl;
					std::exit(1);
				}
			}
			// Other titratable amino acids
			else if (proteinAcceptingMap.find(hydrogen->res_name) != proteinAcceptingMap.end()) {
				parentAtomName = proteinAcceptingMap.at(hydrogen->res_name);
			}
			else {
				std::cerr << "ERROR: Donor H not recognized as titratable residue." << std::endl;
				std::cerr << hydrogen->toString() << std::endl;
				std::exit(1);
			}

			parentAtoms.push_back(parentAtomName);
			std::string key = hydrogen->res_id + parentAtomName;
			parentAtomsMap[key] = hydrogen;
		}

		// Find indices to insert new hydrogens
		std::vector<std::pair<int, std::shared_ptr<Atom>>> insertionIndices;
		for (int i = 0; i < atoms.size(); i++) {

			// Check if accepting residue
			if (std::find(parentResidues.begin(), parentResidues.end(), atoms[i]->res_id) != parentResidues.end()) {

				// Check if corresponding accepting atom in accepting residue
				if (std::find(parentAtoms.begin(), parentAtoms.end(), atoms[i]->atom_name) != parentAtoms.end()) {

					// Find key to map to new hydrogen
					std::string key = atoms[i]->res_id + atoms[i]->atom_name;
					if (parentAtomsMap.find(key) != parentAtomsMap.end()) {
						insertionIndices.push_back(std::pair(i, parentAtomsMap.at(key)));
					}
				}
			}

			// Only need protein atoms
			else if (atoms[i]->res_name == "SOL" || atoms[i]->res_name == "NNN") {
				break;
			}
		}
		if (insertionIndices.size() != newHydrogens.size()) {
			std::cerr << "ERROR: Number of insertion indices does not match number of new hydrogens." << std::endl;
			std::exit(1);
		}

		// Sort insertion indice pairs by index in ascending order
		std::sort(insertionIndices.begin(), insertionIndices.end(), [](auto& left, auto& right) {
			return left.first > right.first;
			});

		// Insert hydrogens at specified indices while maintaing order of vector
		for (const auto& [index, atom] : insertionIndices) {
			if (index <= atoms.size()) { // Ensure the index is valid
				atoms.insert(atoms.begin() + index + 1, atom);
			}
			else {
				std::cerr << "Invalid index for insertion into coordInfo.atoms: " << index << std::endl;
			}
		}
	}

	// Facilate Grotthuss mechanism for H3O+. Given that products = reactants, can just swap the two.
	static void doGrotthussH3O(CoordInfo& coordInfo, const Exchange& exchange) {

		auto hydronium = exchange.donor->parent.lock();
		auto water = exchange.acceptor->parent.lock();

		// Ensures position of exchanged H3O+ H to third, as it will remain in new H3O+ product
		if (exchange.donor->atom_name == "HW1") {
			// Indices to swap
			int idx1 = exchange.donor->atom_num;
			int idx2 = exchange.donor->atom_num + 3;
			swapCoordVel(coordInfo.atoms[idx1], coordInfo.atoms[idx2]);
		}
		else if (exchange.donor->atom_name == "HW2") {
			int idx1 = exchange.donor->atom_num;
			int idx2 = exchange.donor->atom_num + 2;
			swapCoordVel(coordInfo.atoms[idx1], coordInfo.atoms[idx2]);
		}
		else if (exchange.donor->atom_name == "HW3") {
			; // Pass
		}
		else {
			std::cerr << "H3O+ donor atom name not recognized, likely error in coordInfo.atoms ordering." << std::endl;
			exit(1);
		}

		// Swap H3O+ and water coordinates (exlcuding the already swapped H)
		for (int i = 0; i < 4; i++) {
			int idx1 = hydronium->atoms[0].lock()->atom_num + i;
			int idx2 = water->atoms[0].lock()->atom_num + i;
			swapCoordVel(coordInfo.atoms[idx1], coordInfo.atoms[idx2]);
		}

		// Correct geometry of newly formed H3O+
		std::vector<std::shared_ptr<Atom>> hydAtoms;
		for (const auto& atom : hydronium->atoms) {
			hydAtoms.push_back(atom.lock());
		}
		H3O_transform(hydAtoms);
		for (size_t i = 0; i < hydAtoms.size(); i++) {
			hydronium->atoms[i].lock()->coord = hydAtoms[i]->coord;
		}


		// Correct geometry of newly formed  water
		std::vector<std::shared_ptr<Atom>> waterAtoms;
		for (const auto& atom : water->atoms) {
			waterAtoms.push_back(atom.lock());
		}
		water_transform(waterAtoms);
		for (size_t i = 0; i < waterAtoms.size(); i++) {
			water->atoms[i].lock()->coord = waterAtoms[i]->coord;
		}
	}

	// Facilate Grotthuss mechanism for OH-. Similar to above, given that products = reactants, can just swap the two.
	static void doGrotthussOH(CoordInfo& coordInfo, const Exchange& exchange) {

		auto hydroxide = exchange.acceptor->parent.lock();
		auto water = exchange.donor->parent.lock();

		// Swap oxygens
		swapCoordVel(coordInfo.atoms[hydroxide->atoms[0].lock()->atom_num], 
			coordInfo.atoms[water->atoms[0].lock()->atom_num]);

		// Swap non-exchanged hydrogens
		if (exchange.donor->atom_name == "HW2") {
			swapCoordVel(coordInfo.atoms[hydroxide->atoms[1].lock()->atom_num], 
				coordInfo.atoms[water->atoms[2].lock()->atom_num]);
		}
		else if (exchange.donor->atom_name == "HW3") {
			swapCoordVel(coordInfo.atoms[hydroxide->atoms[1].lock()->atom_num], 
				coordInfo.atoms[water->atoms[1].lock()->atom_num]);
		}
		else {
			std::cerr << "H2O donor atom name not recognized, likely error in coordInfo.atoms ordering." << std::endl;
			exit(1);
		}

		// Correct geometry of newly formed water
		std::vector<std::shared_ptr<Atom>> waterAtoms;
		for (const auto& atom : water->atoms) {
			waterAtoms.push_back(atom.lock());
		}
		water_transform(waterAtoms);
		for (size_t i = 0; i < waterAtoms.size(); i++) {
			water->atoms[i].lock()->coord = waterAtoms[i]->coord;
		}
	}

	// Schedule H to be deleted to deprotonate and amino acid, adds proton to the atomsToDelete list
	static void deprotonateAA(CoordInfo& coordInfo, std::vector<int>& atomsToDelete, const Exchange& exchange) {

		// Atom index of the donor hydrogen
		int donorIndex = exchange.donor->atom_num;

		// Add to list of atoms to delete
		atomsToDelete.push_back(donorIndex);

		// Ideally the function would end here, but have to account for multiple protons & atom ordering
		// N-termini
		std::vector<std::string> nterm_H_names = { "H1", "H2", "H3" };
		if ((std::find(nterm_H_names.begin(), nterm_H_names.end(), exchange.donor->atom_name) != nterm_H_names.end())) {
			if (exchange.donor->atom_name == "H1") {
				int idx1 = donorIndex + 1;
				int idx2 = donorIndex + 2;
				coordInfo.atoms[idx1]->atom_name = "H1";
				coordInfo.atoms[idx2]->atom_name = "H2";
			}
			else if (exchange.donor->atom_name == "H2") {
				int idx = donorIndex + 1;
				coordInfo.atoms[idx]->atom_name = "H2";
			}
			else if (exchange.donor->atom_name == "H3") {
				; // Pass
			}
			else {
				std::cerr << "N-term deprotonation state not recognized." << std::endl;
				exit(1);
			}
		}

		// C-termini easy
		else if (exchange.donor->atom_name == "HT2") {
			; // Pass
		}

		// LYS
		else if (exchange.donor->res_name.find("LYS") != std::string::npos) {
			if (exchange.donor->atom_name == "HZ1") {
				int idx1 = donorIndex + 1;
				int idx2 = donorIndex + 2;
				coordInfo.atoms[idx1]->atom_name = "HZ1";
				coordInfo.atoms[idx2]->atom_name = "HZ2";
			}
			else if (exchange.donor->atom_name == "HZ2") {
				int idx = donorIndex + 1;
				coordInfo.atoms[idx]->atom_name = "HZ2";
			}
			else if (exchange.donor->atom_name == "HZ3") {
				; // Pass
			}
			else {
				std::cerr << "LYS deprotonation state not recognized." << std::endl;
				exit(1);
			}
		}

		// ARG is a pain, you have to switch a fair number of residues because only NH1 can be deprotonated
		else if (exchange.donor->res_name.find("ARG") != std::string::npos) {
			if (exchange.donor->atom_name == "HH11") {
				int idx = donorIndex + 1;
				coordInfo.atoms[idx]->atom_name = "HH11";
			}
			else if (exchange.donor->atom_name == "HH21") {
				int idx1 = donorIndex - 4;
				int idx2 = donorIndex - 1;
				coordInfo.atoms[idx1]->atom_name = "NH2";
				coordInfo.atoms[idx2]->atom_name = "NH1";
				swapAtoms(coordInfo.atoms, idx1, idx2);

				idx1 = donorIndex + 1;
				idx2 = donorIndex - 3;
				coordInfo.atoms[idx1]->atom_name = "HH11";
				coordInfo.atoms[idx2]->atom_name = "HH22";
				swapAtoms(coordInfo.atoms, idx1, idx2);

				idx1 = donorIndex - 1;
				idx2 = donorIndex - 2;
				coordInfo.atoms[idx2]->atom_name = "HH21";
				swapAtoms(coordInfo.atoms, idx1, idx2);
			}
			else if (exchange.donor->atom_name == "HH22") {
				int idx1 = donorIndex - 5;
				int idx2 = donorIndex - 2;
				coordInfo.atoms[idx1]->atom_name = "NH2";
				coordInfo.atoms[idx2]->atom_name = "NH1";
				swapAtoms(coordInfo.atoms, idx1, idx2);

				idx1 = donorIndex - 1;
				idx2 = donorIndex - 4;
				coordInfo.atoms[idx1]->atom_name = "HH11";
				coordInfo.atoms[idx2]->atom_name = "HH22";
				swapAtoms(coordInfo.atoms, idx1, idx2);

				idx1 = donorIndex - 2;
				idx2 = donorIndex - 3;
				coordInfo.atoms[idx2]->atom_name = "HH21";
				swapAtoms(coordInfo.atoms, idx1, idx2);
			}
			else if (exchange.donor->atom_name == "HH12") {
				; // Pass
			}
			else {
				std::cerr << "ARG deprotonation state not recognized." << std::endl;
				exit(1);
			}
		}

		// HIS also a pain because atom order between protonated/deprotonated states different
		else if (exchange.donor->res_name.find("HIS") != std::string::npos) {

			int Hindex = exchange.donor->atom_num;
			if (exchange.donor->atom_name == "HD1") {
				swapAtoms(coordInfo.atoms, Hindex - 6, Hindex - 1);
				swapAtoms(coordInfo.atoms, Hindex - 5, Hindex - 4);
				swapAtoms(coordInfo.atoms, Hindex - 4, Hindex + 1);
				swapAtoms(coordInfo.atoms, Hindex - 3, Hindex + 2);
				swapAtoms(coordInfo.atoms, Hindex - 2, Hindex + 2);
				swapAtoms(coordInfo.atoms, Hindex - 1, Hindex + 2);
				swapAtoms(coordInfo.atoms, Hindex + 1, Hindex + 2);
			}
			else if (exchange.donor->atom_name == "HE2") {
				swapAtoms(coordInfo.atoms, Hindex - 4, Hindex + 1);
				swapAtoms(coordInfo.atoms, Hindex - 3, Hindex + 2);
				swapAtoms(coordInfo.atoms, Hindex - 1, Hindex + 3);
				swapAtoms(coordInfo.atoms, Hindex + 1, Hindex + 4);
				swapAtoms(coordInfo.atoms, Hindex + 2, Hindex + 3);
				swapAtoms(coordInfo.atoms, Hindex + 3, Hindex + 4);
			}
		}
		// Other amino acids are easy, just delete the H
		else {
			; // Pass
		}
	}

	// Create new H to protonate an amino acid, adds proton to the hydrogensToAdd list
	static void protonateAA(CoordInfo& coordInfo, std::vector<std::shared_ptr<Atom>>& hydrogensToAdd, const Exchange& exchange,
		const std::unordered_map<std::string, std::string>& proteinHydrogensMap) {

		// Get name of new hydrogen, start by checking if exchange involves termini
		std::string newHydrogenName;
		// N-termini
		if ( exchange.acceptor->atom_name == "N" && 
			(proteinHydrogensMap.find(exchange.acceptor->res_id) != proteinHydrogensMap.end()) ) { 
			
			newHydrogenName = proteinHydrogensMap.at(exchange.acceptor->res_id);
		}
		// C-termini
		else if ( (exchange.acceptor->atom_name == "OT2" || exchange.acceptor->atom_name == "OT1") &&
			(proteinHydrogensMap.find(exchange.acceptor->res_id) != proteinHydrogensMap.end()) ) {

			newHydrogenName = proteinHydrogensMap.at(exchange.acceptor->res_id);

			// Enable protonation of either carboxylate oxygen via swapping
			if (exchange.acceptor->atom_name == "OT1") {
				int idx1 = exchange.acceptor->atom_num;
				int idx2 = exchange.acceptor->atom_num + 1;
				coordInfo.atoms[idx1]->atom_name = "OT2";
				coordInfo.atoms[idx2]->atom_name = "OT1";
				swapAtoms(coordInfo.atoms, exchange.acceptor->atom_num, exchange.acceptor->atom_num + 1);
			}
		}
		// Aspartic acid
		else if (exchange.acceptor->res_name.find("ASP") != std::string::npos) {
			newHydrogenName = proteinHydrogensMap.at(exchange.acceptor->res_name);

			// Enable protonation of either carboxylate oxygen via swapping
			if (exchange.acceptor->atom_name == "OD1") {
				int idx1 = exchange.acceptor->atom_num;
				int idx2 = exchange.acceptor->atom_num + 1;
				coordInfo.atoms[idx1]->atom_name = "OD2";
				coordInfo.atoms[idx2]->atom_name = "OD1";
				swapAtoms(coordInfo.atoms, exchange.acceptor->atom_num, exchange.acceptor->atom_num + 1);
			}
		}
		// Glutamic acid
		else if (exchange.acceptor->res_name.find("GLU") != std::string::npos) {
			newHydrogenName = proteinHydrogensMap.at(exchange.acceptor->res_name);

			// Enable protonation of either carboxylate oxygen via swapping
			if (exchange.acceptor->atom_name == "OE1") {
				int idx1 = exchange.acceptor->atom_num;
				int idx2 = exchange.acceptor->atom_num + 1;
				coordInfo.atoms[idx1]->atom_name = "OE2";
				coordInfo.atoms[idx2]->atom_name = "OE1";
				swapAtoms(coordInfo.atoms, exchange.acceptor->atom_num, exchange.acceptor->atom_num + 1);
			}
		}
		// Histidine special because 2 protonation states (HISE or HISD)
		else if (exchange.acceptor->res_name.find("HIS") != std::string::npos) {

			int hisIndex = exchange.acceptor->atom_num;
			if (exchange.acceptor->atom_name == "ND1") {
				newHydrogenName = "HD1";

				// CHARMM changes the atom ordering of HISD/HISE relative to HISP for god knows what reason
				// Have to correct order here for new protonation state or else MD will fail
				swapAtoms(coordInfo.atoms, hisIndex    , hisIndex + 6);
				swapAtoms(coordInfo.atoms, hisIndex + 1, hisIndex + 7);
				swapAtoms(coordInfo.atoms, hisIndex + 2, hisIndex + 7);
				swapAtoms(coordInfo.atoms, hisIndex + 3, hisIndex + 4);
				swapAtoms(coordInfo.atoms, hisIndex + 4, hisIndex + 5);
				swapAtoms(coordInfo.atoms, hisIndex + 5, hisIndex + 6);
				swapAtoms(coordInfo.atoms, hisIndex + 6, hisIndex + 7);

			}
			else if (exchange.acceptor->atom_name == "NE2") {
				newHydrogenName = "HE2";

				// Also correct atom order here
				swapAtoms(coordInfo.atoms, hisIndex - 5, hisIndex + 1);
				swapAtoms(coordInfo.atoms, hisIndex - 4, hisIndex + 2);
				swapAtoms(coordInfo.atoms, hisIndex - 2, hisIndex);
				swapAtoms(coordInfo.atoms, hisIndex - 1, hisIndex + 1);
				swapAtoms(coordInfo.atoms, hisIndex    , hisIndex + 2);
				swapAtoms(coordInfo.atoms, hisIndex + 1, hisIndex + 2);
			}
			else {
				std::cerr << "HIS protonation state not recognized, likely error in coordInfo.atoms ordering." << std::endl;
				exit(1);
			}
		}
		// Other amino acids
		else if (proteinHydrogensMap.find(exchange.acceptor->res_name) != proteinHydrogensMap.end()) {
			newHydrogenName = proteinHydrogensMap.at(exchange.acceptor->res_name);
		}
		else {
			std::cerr << "ERROR: Acceptor residue not recognized as titratable residue." << std::endl;
			std::exit(1);
		}

		// Find coord of new H by translating along axis to acceptor atom
		std::array<float, 3> newHydrogenCoord = setBondLength(exchange.acceptor->coord, exchange.donor->coord, 0.10f);

		// Creat new H
		auto newHydrogen = std::make_shared<Atom>(
			exchange.acceptor->parent.lock(),
			exchange.acceptor->res_id,
			exchange.acceptor->res_num,
			exchange.acceptor->res_name,
			newHydrogenName,
			-1, // Atom number will be assigned later
			newHydrogenCoord,
			exchange.donor->velocity,
			"H",
			exchange.acceptor->chain
		);
		
		// Add new hydrogen to list
		hydrogensToAdd.push_back(newHydrogen);
	}

	// Create Water from OH- (non-Grotthuss) 
	static void createWaterFromOH(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount,
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old hydroxide
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.acceptor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new H3O+ residue
		std::vector<std::shared_ptr<Atom>> newRes;
		newRes.reserve(resAtoms.size() + 2); // 4 water atoms (including virtual site) + 1 new H
		resCount++;
		std::string newResName = "SOL";
		std::string newResID = std::to_string(resCount) + newResName;
		std::string newChain = resAtoms[0].lock()->chain;	

		std::vector<std::string> newAtomNames = { "OW1", "HW2" };
		for (int i = 0; i < resAtoms.size(); i++) {
			auto resAtom = resAtoms[i].lock();
			auto newAtom = std::make_shared<Atom>(nullptr, newResID, resCount, newResName, 
				newAtomNames[i], 0, resAtom->coord, resAtom->velocity, resAtom->element, resAtom->chain);
			newRes.push_back(newAtom);
		}

		// Add in new Hydrogen
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "HW3", 0,
			setBondLength(exchange.acceptor->coord, exchange.donor->coord, 0.09686f), 
			exchange.donor->velocity, "H", newChain)
		);

		// Add in new virtual site, dummy coord/vel for now, will be corrected during transform
		std::array<float, 3 > dummy = { 0.f, 0.f, 0.f };
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "MW4", 0, dummy, dummy, "M", newChain)
		);

		// Ensure proper geometry of newly formed water
		water_transform(newRes);

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Create H3O+ from water (non-Grotthuss) 
	static void createH3O(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount, 
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old water
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.acceptor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new H3O+ residue
		std::vector<std::shared_ptr<Atom>> newRes;
		newRes.reserve(resAtoms.size() + 1); // 4 water atoms + 1 new H
		resCount++;
		std::string newResName = "HHO";
		std::string newResID = std::to_string(resCount) + newResName;

		std::vector<std::string> newAtomNames = { "OW", "HW1", "HW2", "MW" };
		for (int i = 0; i < resAtoms.size(); i++) {
			auto resAtom = resAtoms[i].lock();
			auto newAtom = std::make_shared<Atom>(nullptr, newResID, resCount, newResName,
				newAtomNames[i], 0, resAtom->coord, resAtom->velocity, resAtom->element, resAtom->chain);
			newRes.push_back(newAtom);
		}

		// Add in new Hydrogen
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "HW3", 0,
			setBondLength(exchange.acceptor->coord, exchange.donor->coord, 0.09686f), 
			exchange.donor->velocity, "H", "Z")
		);

		// Ensure proper geometry of newly formed H3O+
		H3O_transform(newRes);

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Create acetic acid from acetate
	static void createAcetic(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount,
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old 
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.acceptor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new 
		std::vector<std::shared_ptr<Atom>> newRes;
		resCount++;
		std::string newResName = "AHX";
		std::string newResID = std::to_string(resCount) + newResName;
		std::string newChain = resAtoms[0].lock()->chain;
		newRes.reserve(resAtoms.size() + 1);

		std::vector<std::string> newAtomNames = { "C2", "C1", "H21", "H22", "H23", "O2", "O1" };
		for (int i = 0; i < resAtoms.size(); i++) {
			auto resAtom = resAtoms[i].lock();
			auto newAtom = std::make_shared<Atom>(nullptr, newResID, resCount, newResName,
				newAtomNames[i], 0, resAtom->coord, resAtom->velocity, resAtom->element, resAtom->chain);
			newRes.push_back(newAtom);
		}

		// Protonation of either carboxylate O
		if (exchange.acceptor->atom_name == "O1") {
			newRes[5]->atom_name = "O1";
			newRes[6]->atom_name = "O2";
			std::swap(newRes[5], newRes[6]);
		}

		// Add in new Hydrogen
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "HO1", 0,
			setBondLength(exchange.acceptor->coord, exchange.donor->coord, 0.10f), 
			exchange.donor->velocity, "H", newChain)
		);

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Create ammonium from ammonia
	static void createAmmonium(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount,
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.acceptor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new
		std::vector<std::shared_ptr<Atom>> newRes;
		resCount++;
		std::string newResName = "NXH";
		std::string newResID = std::to_string(resCount) + newResName;
		newRes.reserve(resAtoms.size() + 1);

		// Add in new Hydrogen
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "HZ1", 0,
			setBondLength(exchange.acceptor->coord, exchange.donor->coord, 0.10f), 
			exchange.donor->velocity, "H", exchange.acceptor->chain)
		);

		// Add in other atoms
		std::vector<std::string> newAtomNames = { "NZ", "HZ2", "HZ3", "HZ4" };
		for (int i = 0; i < resAtoms.size(); i++) {
			auto resAtom = resAtoms[i].lock();
			auto newAtom = std::make_shared<Atom>(nullptr, newResID, resCount, newResName,
				newAtomNames[i], 0, resAtom->coord, resAtom->velocity, resAtom->element, resAtom->chain);
			newRes.push_back(newAtom);
		}

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Create Water from H3O+ (non-Grotthuss) 
	static void createWaterFromH3O(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount,
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old 
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.donor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new 
		std::vector<std::shared_ptr<Atom>> newRes;
		newRes.reserve(resAtoms.size() - 1); // Subtract donated H
		resCount++;
		std::string newResName = "SOL";
		std::string newResID = std::to_string(resCount) + newResName;

		// Add in new O
		auto oxygen = resAtoms[0].lock();
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "OW1", 0, oxygen->coord, 
			oxygen->velocity , oxygen->element, oxygen->chain)
		);

		// Add in hydrogens, order depends on which H being donated
		std::shared_ptr<Atom> firstH;
		std::shared_ptr<Atom> secondH;
		if (exchange.donor->atom_name == "HW1") {
			firstH = resAtoms[2].lock();
			secondH = resAtoms[4].lock();
		}
		else if (exchange.donor->atom_name == "HW2") {
			firstH = resAtoms[1].lock();
			secondH = resAtoms[4].lock();
		}
		else if (exchange.donor->atom_name == "HW3") {
			firstH = resAtoms[1].lock();
			secondH = resAtoms[2].lock();
		}
		else {
			std::cerr << "H3O+ non-donor H name not recognized, likely error in coordInfo.atoms ordering." << std::endl;
			exit(1);
		}

		// First hydrogen
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "HW2", 0, firstH->coord, 
			firstH->velocity, "H", firstH->chain)
		);

		// Second hydrogen
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "HW3", 0, secondH->coord,
			secondH->velocity, "H", secondH->chain)
		);

		// Virtual site
		auto vsite = resAtoms[3].lock();
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "MW4", 0, vsite->coord,
			vsite->velocity, "M", vsite->chain)
		);

		// Ensure proper geometry of newly formed water
		water_transform(newRes);

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Create OH- from water (non-Grotthuss) 
	static void createOH(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount,
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old 
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.donor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new 
		std::vector<std::shared_ptr<Atom>> newRes;
		newRes.reserve(2); // Subtract donated H and virtual site
		resCount++;
		std::string newResName = "OHX";
		std::string newResID = std::to_string(resCount) + newResName;

		// Add in new O
		auto oxygen = resAtoms[0].lock();
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "O1", 0, oxygen->coord, 
			oxygen->velocity, oxygen->element, oxygen->chain)
		);

		// Add in new H, take care which one
		int Hindex;
		if (exchange.donor->atom_name == "HW2") {
			Hindex = 2;
		}
		else if (exchange.donor->atom_name == "HW3") {
			Hindex = 1;
		}
		else {
			std::cerr << "H2O donor atom name not recognized, likely error in coordInfo.atoms ordering." << std::endl;
			exit(1);
		}

		auto hydrogen = resAtoms[Hindex].lock();
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "H1", 0, hydrogen->coord, hydrogen->velocity, 
			hydrogen->element, hydrogen->chain)
		);

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Create acetate from acetic acid
	static void createAcetate(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount,
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old 
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.donor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new 
		std::vector<std::shared_ptr<Atom>> newRes;
		newRes.reserve(resAtoms.size() - 1); // Subtract donated H
		resCount++;
		std::string newResName = "ATX";
		std::string newResID = std::to_string(resCount) + newResName;

		std::vector<std::string> newAtomNames = { "C1", "C2", "H1", "H2", "H3", "O1", "O2" };
		for (int i = 0; i < resAtoms.size() - 1; i++) {
			auto resAtom = resAtoms[i].lock();
			auto newAtom = std::make_shared<Atom>(nullptr, newResID, resCount, newResName,
				newAtomNames[i], 0, resAtom->coord, resAtom->velocity, resAtom->element, resAtom->chain);
			newRes.push_back(newAtom);
		}

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Create ammonia from ammonium
	static void createAmmonia(std::vector<std::shared_ptr<Atom>>& productsToAdd, std::vector<int>& atomsToDelete, int& resCount,
		const CoordInfo& coordInfo, const Exchange& exchange) {

		// Delete old 
		std::vector<std::weak_ptr<Atom>> resAtoms = exchange.donor->parent.lock()->atoms;
		for (const auto& atom : resAtoms) {
			atomsToDelete.push_back(atom.lock()->atom_num);
		}

		// Begin creating new 
		std::vector<std::shared_ptr<Atom>> newRes;
		newRes.reserve(resAtoms.size() - 1); // Subtract donated H
		resCount++;
		std::string newResName = "NXX";
		std::string newResID = std::to_string(resCount) + newResName;

		// Add in new N
		auto nitrogen = resAtoms[1].lock();
		newRes.push_back(std::make_shared<Atom>(
			nullptr, newResID, resCount, newResName, "N1", 0, nitrogen->coord, 
			nitrogen->velocity, nitrogen->element, nitrogen->chain)
		);

		// Find which H's to keep 
		std::vector<std::shared_ptr<Atom>> nonDonatedHs;
		if (exchange.donor->atom_name == "HZ1") {
			nonDonatedHs.push_back(resAtoms[2].lock());
			nonDonatedHs.push_back(resAtoms[3].lock());
			nonDonatedHs.push_back(resAtoms[4].lock());
		}
		else if (exchange.donor->atom_name == "HZ2") {
			nonDonatedHs.push_back(resAtoms[0].lock());
			nonDonatedHs.push_back(resAtoms[3].lock());
			nonDonatedHs.push_back(resAtoms[4].lock());
		}
		else if (exchange.donor->atom_name == "HZ3") {
			nonDonatedHs.push_back(resAtoms[0].lock());
			nonDonatedHs.push_back(resAtoms[2].lock());
			nonDonatedHs.push_back(resAtoms[4].lock());
		}
		else if (exchange.donor->atom_name == "HZ4") {
			nonDonatedHs.push_back(resAtoms[0].lock());
			nonDonatedHs.push_back(resAtoms[2].lock());
			nonDonatedHs.push_back(resAtoms[3].lock());
		}
		else {
			std::cerr << "NH3 donor atom name not recognized, likely error in coordInfo.atoms ordering." << std::endl;
			exit(1);
		}

		std::vector<std::string> newAtomNames = { "H11", "H12", "H13" };
		for (int i = 0; i < 3; i++) {
			std::shared_ptr<Atom> resAtom = nonDonatedHs[i];
			auto newAtom = std::make_shared<Atom>(nullptr, newResID, resCount, newResName,
				newAtomNames[i], 0, resAtom->coord, resAtom->velocity, resAtom->element, resAtom->chain);
			newRes.push_back(newAtom);
		}

		// Schedule for insertion by adding to products
		productsToAdd.insert(productsToAdd.end(), newRes.begin(), newRes.end());
	}

	// Final facilitation of all exchanges by coordinating the above functions. Modifies coordInfo object in place.
	void doExchanges(CoordInfo& coordInfo, const std::vector<Exchange>& exchanges, const bool& pairs, const bool& prot,
		const std::vector<std::string>& topOrder, const Config& config) {

		if (!pairs) {
			return;
		}

		// Lists associated with changing protonation states
		int resCount = coordInfo.residues.size(); // Running tally of residues
		std::vector<std::shared_ptr<Atom>> hydrogensToAdd; // New protein H's
		std::vector<int> atomsToDelete; // Old reactant residues deleted
		std::vector<std::shared_ptr<Atom>> productsToAdd; // New product residues added

		// Because protein exchanges handled differently, store useful protein info. First, get protein residue names, including termini
		std::vector<std::string> proteinResNames = getProteinResNames(coordInfo);

		// If protonating an amino acid, name of the new H per each amino acid residue
		// HIS, NTERM, and CTERM are special cases, so handled seperately
		std::unordered_map<std::string, std::string> proteinHydrogensMap = {
			{ "LYS", "HZ3" }, { "ARG", "HH12" }, { "ASP", "HD2" }, { "GLU", "HE2" }
		};

		// This stores the atom name in front of the H names above for insertion to Atom list
		std::unordered_map<std::string, std::string> proteinAcceptingMap = {
			{ "LYS", "HZ2" }, { "ARG", "HH11" }, { "ASP", "OD2" }, { "GLU", "OE2" },
		};
		
		// Termini info
		std::vector<std::string> ntermResIDs;
		std::vector<std::string> ctermResIDs;
		for (const auto& monomer : coordInfo.proteins) {

			// Store N and C termini resID's
			std::string ntermResID = monomer->residues[0].lock()->res_id;
			std::string ctermResID = monomer->residues.back().lock()->res_id;
			ntermResIDs.push_back(ntermResID);
			ctermResIDs.push_back(ctermResID);

			// Add termini to proteinHydrogensMap
			// Using resID allows differentation of termini if protein a complex
			proteinHydrogensMap[ntermResID] = "H3";
			proteinHydrogensMap[ctermResID] = "HT2";

			// Add termini to proteinAcceptingMap
			proteinAcceptingMap[ntermResID] = "H2";
			proteinAcceptingMap[ctermResID] = "OT2";
		}

		// Change donor protonation states
		for (const auto& exchange : exchanges) {

			// Grotthuss mechanism for H3O+
			if (exchange.donor->res_name == "HHO" && exchange.acceptor->res_name == "SOL") {
				doGrotthussH3O(coordInfo, exchange);
			}

			// Grotthuss mechanism for OH-
			else if (exchange.acceptor->res_name == "OHX" && exchange.donor->res_name == "SOL") {
				doGrotthussOH(coordInfo, exchange);
			}

			// Deprotonated water and form OH- (non-Grotthuss)
			else if (exchange.donor->res_name == "SOL") {
				createOH(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Deprotonate H3O+ and form water (non-Grotthuss)
			else if (exchange.donor->res_name == "HHO") {
				createWaterFromH3O(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Deprotonated acetic acid
			else if (exchange.donor->res_name == "AHX") {
				createAcetate(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Deprotonate ammonium
			else if (exchange.donor->res_name == "NXH") {
				createAmmonia(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Deprotonate amino acid
			else if (std::find(proteinResNames.begin(), proteinResNames.end(), exchange.donor->res_name) != proteinResNames.end()) {
				deprotonateAA(coordInfo, atomsToDelete, exchange);
			}

			else {
				std::cerr << "Donor residue name not recognized." << std::endl;
			}
		}

		// Change acceptor protonation states
		for (const auto& exchange : exchanges) {

			// Grotthuss mechanism are special, handled in donor loop
			if ( (exchange.donor->res_name == "HHO" && exchange.acceptor->res_name == "SOL") || 
				 (exchange.acceptor->res_name == "OHX" && exchange.donor->res_name == "SOL")) {
				continue;
			}

			// Form H3O+ from water (non-Grotthuss)
			else if (exchange.acceptor->res_name == "SOL") {
				createH3O(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Form water from OH- (non-Grotthuss)
			else if (exchange.acceptor->res_name == "OHX") {
				createWaterFromOH(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Protonate acetate
			else if (exchange.acceptor->res_name == "ATX") {
				createAcetic(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Protonate ammonia
			else if (exchange.acceptor->res_name == "NXX") {
				createAmmonium(productsToAdd, atomsToDelete, resCount, coordInfo, exchange);
			}

			// Protonate amino acid
			else if (std::find(proteinResNames.begin(), proteinResNames.end(), exchange.acceptor->res_name) != proteinResNames.end()) {
				protonateAA(coordInfo, hydrogensToAdd, exchange, proteinHydrogensMap);
			}

			else {
				std::cerr << "ERROR: Acceptor residue name not recognized." << std::endl;
				std::cerr << exchange.toString() << std::endl;
				exit(1);
			}
		}

		// Begin updating coordinate/topology files, start by deleting old reactants
		deleteAtoms(coordInfo.atoms, atomsToDelete);

		// If protonating any amino acids, add new H's here
		if (prot) {
			insertAtoms(coordInfo.atoms, hydrogensToAdd, proteinAcceptingMap);
		}
		
		// Add new products
		for (const auto& atom : productsToAdd) {
			coordInfo.atoms.push_back(atom);
		}

		// Update coordInfo, including proper ordering for .top
		reorderAtoms(coordInfo, topOrder);
	}

}