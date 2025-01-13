#include "Core.h"
#include <random>
#include <set>
#include <algorithm>
#include <functional>

namespace Core {
	
	// Get residue names of all amino acids, including termini.
	std::vector<std::string> getProteinResNames(const CoordInfo& coordInfo) {
		std::vector<std::string> proteinResNames = { "ARG", "ASP", "GLU", "HIS", "LYS" };
		for (const auto& monomer : coordInfo.proteins) {
			auto ntermResidue = monomer->residues.front().lock();
			auto ctermResidue = monomer->residues.back().lock();

			if (ntermResidue) {
				std::string ntermResName = ntermResidue->res_name;
				if (std::find(proteinResNames.begin(), proteinResNames.end(), ntermResName) == proteinResNames.end()) {
					proteinResNames.push_back(ntermResName);
				}
			}

			if (ctermResidue) {
				std::string ctermResName = ctermResidue->res_name;
				if (std::find(proteinResNames.begin(), proteinResNames.end(), ctermResName) == proteinResNames.end()) {
					proteinResNames.push_back(ctermResName);
				}
			}
		}

		return proteinResNames;
	}

	// Gas phase correction that smoothly switches from gas phase basicity to pka based on number of waters
	// Scaling factor from A. Kumar et al. Physical Chemistry Chemical Physics 2022 Vol. 24 Issue 30 Pages 18236-18244
	static float gasPhaseCorrection(const TitratableAtom& atom, const float& temperature, const int& nearWaters) {
		float solE = 0.019144 * atom.pka * temperature; // Energy from pKa
		return ((atom.gpb - solE) * std::exp(-0.30312 * nearWaters)) + solE;
	}

	// Computes energy of exchange considering degree of solvation
	static float computeEnergy(const TitratableAtom& donor, const TitratableAtom& acceptor, const float& temperature, 
		const int& nearWaters) {
		float energy;

		// If waters > 30, use pKa only
		if (nearWaters > 30) {
			energy = 0.019144 * (donor.pka - acceptor.pka) * temperature;
		}	

		// Gas phase correction
		else {
			float donorE = gasPhaseCorrection(donor, temperature, nearWaters);
			float acceptorE = gasPhaseCorrection(acceptor, temperature, nearWaters);
			energy =  donorE - acceptorE;
		}

		return energy;
	}

	// Metropolis Criterion Monte Carlo (MCMC) sampling
	static bool MCMC(const float& deltaE, const float& temperature) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> dis(0.0, 1.0);
		float randNum = dis(gen);
		float boltzmann = std::exp(-deltaE / (0.00831446261815324 * temperature));
		return randNum < boltzmann;
	}

	// Determines if a potential exchange is a Grotthuss exchange
	static bool isGrotthuss(const Exchange& exchange) {
		return (exchange.donor->res_name == "HHO" && exchange.acceptor->res_name == "SOL") ||
			(exchange.acceptor->res_name == "OHX" && exchange.donor->res_name == "SOL");
	}

	/* Determines whether exchange changes pH of solvating cluster, if so, only one pH changing 
		exchange allowed, per cluster, per step. Grotthuss exchanges do not change pH. */
	static bool isAllowedExchange(const Exchange& exchange) {
		bool isAllowed = true;

		std::vector<std::shared_ptr<Atom>> exchangeAtoms = { exchange.donor, exchange.acceptor };
		for (const auto& atom : exchangeAtoms) {
			if ( (atom->res_name == "HHO" || atom->res_name == "OHX" || atom->res_name == "SOL") && !(isGrotthuss(exchange)) ) {
				isAllowed = false;
			}
		}

		return isAllowed;
	}

	// Delets potential exchanges from vector based on inputted indices
	static void deleteExchanges(std::vector<Exchange>& exchanges, const std::set<int>& indices) {
		// Iterate through the set in reverse order to avoid invalidating indices
		for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
			if (*it >= 0 && *it < static_cast<int>(exchanges.size())) {
				exchanges.erase(exchanges.begin() + *it);
			}
		}
	}

	// Find potential exchange donor/acceptor pairs where one residue is water
	static void findWaterPairs(std::vector<Exchange>& exchanges, const CoordInfo& coordInfo, const std::vector<TitratableAtom>& nonWaterAtoms,
		const std::vector<int>& clusterIDs, const std::unordered_map<int, Cluster>& clusters, const float& temperature,
		const std::vector<float>& charges, const int& step, const int& hop, const bool& isDonor) {

		for (const auto& atom : nonWaterAtoms) {

			// Non-Grotthuss exchanges involving water only on hop 0
			if (atom.atom->res_name == "HHO" || atom.atom->res_name == "OHX") {
				;
			}
			else if (hop > 0) {
				continue;
			}

			// Find closest water oxygen (if donor) or hydrogen (if acceptor)
			int minIndex;
			float minSquaredDist = 9999.9f;
			float minDist = 9999.9f;
			if (isDonor) {
				for (int i = 0; i < coordInfo.waterOCoords.size(); i++) {
					const float dx = coordInfo.waterOCoords[i][0] - atom.coord[0];
					const float dy = coordInfo.waterOCoords[i][1] - atom.coord[1];
					const float dz = coordInfo.waterOCoords[i][2] - atom.coord[2];
					const float distanceSquared = dx * dx + dy * dy + dz * dz;

					if (distanceSquared < minSquaredDist) {
						minIndex = i;
						minSquaredDist = distanceSquared;
						minDist = std::sqrt(distanceSquared);
					}
				}
			}
			else {
				for (int i = 0; i < coordInfo.waterHCoords.size(); i++) {
					const float dx = coordInfo.waterHCoords[i][0] - atom.coord[0];
					const float dy = coordInfo.waterHCoords[i][1] - atom.coord[1];
					const float dz = coordInfo.waterHCoords[i][2] - atom.coord[2];
					const float distanceSquared = dx * dx + dy * dy + dz * dz;

					if (distanceSquared < minSquaredDist) {
						minIndex = i;
						minSquaredDist = distanceSquared;
						minDist = std::sqrt(distanceSquared);
					}
				}
			}

			// Only consider exchanges within 0.25nm of each other
			if (minDist > 0.25f) {
				continue;
			}

			// Compute energy and add to potential exchanges 
			else {

				// For formation of H3O+
				if (isDonor) {
					// Get atom object from min index and pack into TitratableAtom object for E calculation
					std::shared_ptr<Atom> waterAtom = coordInfo.residueMap.at("SOL")[minIndex]->atoms[0].lock();
					int clusterID = clusterIDs[minIndex];
					Cluster waterCluster = clusters.at(clusterID);

					// Set pKa to pH of the cluster, and set GPB equal to that of H3O+
					TitratableAtom acceptor = { waterAtom, waterAtom->coord, waterCluster.numWaters, clusterID, waterCluster.pH, 659.82f };

					// Energies for Grotthuss exchanges are computed differently
					float energy;
					int nearWaters;
					if (atom.nearWaters >= acceptor.nearWaters) {
						nearWaters = atom.nearWaters;
					}
					else {
						nearWaters = acceptor.nearWaters;
					}
					if (atom.atom->res_name == "HHO") {
						energy = grotthussEnergy(coordInfo, charges, atom, acceptor, hop);
					}
					else {
						energy = computeEnergy(atom, acceptor, temperature, nearWaters);
					}

					exchanges.push_back(Exchange{ atom.atom, acceptor.atom, energy, step, hop, atom.clusterID, waterCluster.pH, nearWaters });
				}

				// For formation of OH-
				else {
					// Get atom object from min index and pack into TitratableAtom object for E calculation
					int waterIndex = minIndex / 2; // Account for 2 hydrogens per water
					int Hindex = (minIndex % 2) + 1; // For indexing within water residue
					std::shared_ptr<Atom> waterAtom = coordInfo.residueMap.at("SOL")[waterIndex]->atoms[Hindex].lock();
					int clusterID = clusterIDs[waterIndex];
					Cluster waterCluster = clusters.at(clusterID);

					// Set pKa to pH of the cluster, and set GPB equal to that of OH-
					TitratableAtom donor = { waterAtom, waterAtom->coord, waterCluster.numWaters, clusterID, waterCluster.pH, 1605.40f };

					// Energies for Grotthuss exchanges are computed differently
					float energy;
					int nearWaters;
					if (atom.nearWaters >= donor.nearWaters) {
						nearWaters = atom.nearWaters;
					}
					else {
						nearWaters = donor.nearWaters;
					}
					if (atom.atom->res_name == "OHX") {
						energy = grotthussEnergy(coordInfo, charges, donor, atom, hop);
					}
					else {
						energy = computeEnergy(donor, atom, temperature, nearWaters);
					}
					exchanges.push_back(Exchange{ donor.atom, atom.atom, energy, step, hop, atom.clusterID, waterCluster.pH, nearWaters });
				}
			}
		}
	}

	// Find exchange donor/acceptor pairs to facilitate via energy calculations and pair selection criteria
	std::vector<Exchange> findPairs(const CoordInfo& coordInfo, const TitratableSites& titSites,
		const std::vector<int>& clusterIDs, const std::unordered_map<int, Cluster>& clusters, const float& temperature, 
		const int& step, const int& hop, const std::vector<float>& charges, std::vector<std::array <float, 3>>& skipCoords,
		bool& pairs, bool& prot) {

		// Get names of all amino acids, including termini
		std::vector<std::string> proteinResNames = getProteinResNames(coordInfo);

		// Start with exchanges not involving water first
		std::vector<Exchange> exchanges;
		for (const auto& donor : titSites.donors) {
			for (const auto& acceptor : titSites.acceptors) {

				// No self-exchanges 
				if (donor.atom->res_num == acceptor.atom->res_num) {
					continue;
				}

				// No intramolecular protein exchanges, ie, mobile protons (for now!)
				if (std::find(proteinResNames.begin(), proteinResNames.end(), donor.atom->res_name) != proteinResNames.end()
					&& std::find(proteinResNames.begin(), proteinResNames.end(), acceptor.atom->res_name) != proteinResNames.end()) {
					continue;
				}

				// Only consider exchanges within 0.25nm of each other
				float dist = getDistance(donor.atom->coord, acceptor.atom->coord);
				if (dist > 0.25f) {
					continue;
				}

				// Compute energy and add to potential exchanges
				else {
					int nearWaters;
					if (donor.nearWaters >= acceptor.nearWaters) {
						nearWaters = donor.nearWaters;
					}
					else {
						nearWaters = acceptor.nearWaters;
					}
					float energy = computeEnergy(donor, acceptor, temperature, nearWaters);
					Cluster cluster = clusters.at(donor.clusterID);
					exchanges.push_back(Exchange{ donor.atom, acceptor.atom, energy, step, hop, donor.clusterID, cluster.pH, nearWaters });
				}
			}
		}

		// Exchanges involving water molecule
		if (coordInfo.numResidues.at("SOL") > 0) {
			findWaterPairs(exchanges, coordInfo, titSites.donors, clusterIDs, clusters, temperature, 
				charges, step, hop, true); // Water as accepting site
			findWaterPairs(exchanges, coordInfo, titSites.acceptors, clusterIDs, clusters, temperature, 
				charges, step, hop, false); // Water as H donor
		}

		// Begin selecting final, accepted donor/acceptor pairs based on multiple criteria, start with MCMC sampling
		std::set<int> pairsToDelete;
		for (int i = 0; i < exchanges.size(); i++) {
			if (!MCMC(exchanges[i].energy, temperature)) {
				pairsToDelete.insert(i);
			}
		}
		deleteExchanges(exchanges, pairsToDelete);

		// Ensure no two exchanges have exact identical energies to enable comparison
		size_t numExchanges = exchanges.size();
		if (numExchanges > 1) {
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<float> dis(-0.01, 0.01);

			std::set<float> uniqueNumbers;
			while (uniqueNumbers.size() < numExchanges) {
				uniqueNumbers.insert(dis(gen));
			}
			std::vector<float> result(uniqueNumbers.begin(), uniqueNumbers.end());

			for (size_t i = 1; i < numExchanges; i++) {
				exchanges[i].energy += result[i];
			}
		}

		std::vector<float> energies;
		for (auto& exchange : exchanges) {
			if (std::find(energies.begin(), energies.end(), exchange.energy) != energies.end()) {
				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_real_distribution<float> dis(-0.01, 0.01);
				exchange.energy += dis(gen);
			}
			energies.push_back(exchange.energy);
		}

		// Cannot protonate a carboxylate/Arg if angle ~180 degrees or else simulation breaks
		pairsToDelete.clear();
		for (int i = 0; i < exchanges.size(); i++) {
			
			// Find name of carbon carboxylate if acceptor is a carboxylate group
			std::string carbonName;
			if (exchanges[i].acceptor->res_name == "ATX") {
				carbonName = "C2";
			} 
			else if (exchanges[i].acceptor->res_name == "GLU") {
				carbonName = "CD";
			}
			else if (exchanges[i].acceptor->res_name == "ASP") {
				carbonName = "CG";
			}
			else if (exchanges[i].acceptor->res_name == "ARG") {
				carbonName = "CZ";
			}
			else {
				continue;
			}

			// Find carbon atom in carboxylate group and compute angle
			for (const auto& weak_atom : exchanges[i].acceptor->parent.lock()->atoms) {
				auto atom = weak_atom.lock();
				// FIX - What should 180* bug cutoff be if anything?
				if (atom->atom_name == carbonName) {
					if (std::abs(getAngle(atom->coord, exchanges[i].acceptor->coord, exchanges[i].donor->coord)) < 3.0f) {
						pairsToDelete.insert(i);
					}
					break;
				}
			}
		}
		deleteExchanges(exchanges, pairsToDelete);

		/* Prevent 'rattling' ie, the passing of proton of between two residues over multiple hops.
			Because atoms are updated after each hop, can't identifiers like res_id, must use coordinates. */
		pairsToDelete.clear();
		if (skipCoords.size() > 0) {
			for (int i = 0; i < exchanges.size(); i++) {
				std::vector<float> distances = norm(skipCoords, exchanges[i].acceptor->coord);
				if (*std::min_element(distances.begin(), distances.end()) < 0.10f) {
					pairsToDelete.insert(i);
				}
			}
			deleteExchanges(exchanges, pairsToDelete);
		}

		// Find most energitically favorable Grotthuss exchange per H3O+ or OH-, water pair
		std::unordered_map<std::string, float> grotthussEnergies;
		for (const auto& exchange : exchanges) {
			if (isGrotthuss(exchange)) {
				std::vector<std::shared_ptr<Atom>> exchangeAtoms = { exchange.donor, exchange.acceptor };
				for (const auto& atom : exchangeAtoms) {
					if (grotthussEnergies.find(atom->res_id) == grotthussEnergies.end()) {
						grotthussEnergies[atom->res_id] = exchange.energy;
					}
					else {
						if (exchange.energy < grotthussEnergies[atom->res_id]) {
							grotthussEnergies[atom->res_id] = exchange.energy;
						}
					}
				}
			}
		}

		// Remove less favorable Grotthuss residue exchanges
		pairsToDelete.clear();
		for (int i = 0; i < exchanges.size(); i++) {
			if (isGrotthuss(exchanges[i])) {
				if (exchanges[i].energy > grotthussEnergies[exchanges[i].donor->res_id]
					|| exchanges[i].energy > grotthussEnergies[exchanges[i].acceptor->res_id]) {
					pairsToDelete.insert(i);
				}
			}
		}
		deleteExchanges(exchanges, pairsToDelete);

		// Find most energitically favorable exchange per residue 
		std::unordered_map<std::string, float> exchangeEnergies;
		for (auto& exchange : exchanges) {

			/* Set Grotthuss exchanges to have an energy to zero for comparison purposes
				with other non-Grotthuss exchanges as the delta G is zero as products = reactants */
			if (isGrotthuss(exchange)) {
				// Rand number to break ties
				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_real_distribution<float> dis(0.0, 1.0);
				float randNum = dis(gen);
				exchange.energy = 0.0f + randNum;
			}

			std::vector<std::shared_ptr<Atom>> exchangeAtoms = { exchange.donor, exchange.acceptor };
			for (const auto& atom : exchangeAtoms) {
				if (exchangeEnergies.find(atom->res_id) == exchangeEnergies.end()) {
					exchangeEnergies[atom->res_id] = exchange.energy;
				}
				else {
					if (exchange.energy < exchangeEnergies[atom->res_id]) {
						exchangeEnergies[atom->res_id] = exchange.energy;
					}
				}
			}
		}

		// Remove less favorable residue exchanges
		pairsToDelete.clear();
		for (int i = 0; i < exchanges.size(); i++) {
			if (exchanges[i].energy > exchangeEnergies[exchanges[i].donor->res_id] 
				|| exchanges[i].energy > exchangeEnergies[exchanges[i].acceptor->res_id]) {
				pairsToDelete.insert(i);
			}
		}
		deleteExchanges(exchanges, pairsToDelete);

		// Find most energitically favorable pH changing exchange, per water cluster
		std::unordered_map<int, float> clusterEnergies;
		for (auto& exchange : exchanges) {

			// Determine if pH-changing exchange
			if (!isAllowedExchange(exchange)) {

				if (clusterEnergies.find(exchange.clusterID) == clusterEnergies.end()) {
					clusterEnergies[exchange.clusterID] = exchange.energy;
				}
				else {
					if (exchange.energy < clusterEnergies[exchange.clusterID]) {
						clusterEnergies[exchange.clusterID] = exchange.energy;
					}
				}
			}
		}

		// Remove less favorable pH changing exchanges
		pairsToDelete.clear();
		for (int i = 0; i < exchanges.size(); i++) {
			if (!isAllowedExchange(exchanges[i])) {
				if (exchanges[i].energy > clusterEnergies[exchanges[i].clusterID]) {
					pairsToDelete.insert(i);
				}
			}
		}
		deleteExchanges(exchanges, pairsToDelete);

		/* Update skipCoords to prevent 'rattling' of protons between donor/acceptor pair as well as 
			flags that determine if any pairs accepted(pairs bool) or if new protein top needed (prot bool) */
		if (exchanges.size() > 0) {
			pairs = true;
			for (const auto& exchange : exchanges) {
				std::vector<std::shared_ptr<Atom>> exchangeAtoms = { exchange.donor, exchange.acceptor };
				for (const auto& atom : exchangeAtoms) {
					if (std::find(proteinResNames.begin(), proteinResNames.end(), atom->res_name) != proteinResNames.end()) {
						prot = true;
					}
					skipCoords.push_back(atom->coord);
				}
			}
		}

		return exchanges;
	}
	
}