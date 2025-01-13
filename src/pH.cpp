#include "Core.h"
#include <unordered_map>
#include <algorithm>

namespace Core {

	// Find number of waters in each unique cluster
	std::unordered_map<int, int> findClusterWaters(const std::vector<int>& clusterIDs) {
		std::unordered_map<int, int> clusterWaters;
		for (const int& cluster : clusterIDs) {
			clusterWaters[cluster]++;
		}
		return clusterWaters;
	}

	// Via distance, pin residues to the cluster of the closest water
	std::vector<int> pinCluster(const std::vector<std::shared_ptr<Atom>>& titAtoms, const std::vector<std::array<float, 3>>& waterCoords,
		const std::vector<int>& clusterIDs, const std::unordered_map<int, int>& clusterWaters, const bool& donor) {
		std::vector<int> pinnedClusters;

		// If everything in gas phase
		if (clusterIDs.size() == 0) {
			for (int i = 0; i < titAtoms.size(); i++) {
				pinnedClusters.push_back(-1);
			}
			return pinnedClusters;
		}

		// Apply distance cutoff of 0.35 nm to each atom in residue
		for (const auto& titAtom : titAtoms) {

			// Find all waters within 1.00 nm of first titAtom
			std::vector<int> closeWater_indices;
			std::vector<std::array<float, 3>> closeWater_coords;
			for (int i = 0; i < waterCoords.size(); i++) {
				const float dx = waterCoords[i][0] - titAtom->coord[0];
				const float dy = waterCoords[i][1] - titAtom->coord[1];
				const float dz = waterCoords[i][2] - titAtom->coord[2];
				const float distanceSquared = dx * dx + dy * dy + dz * dz;

				if (distanceSquared < 1.00) {
					closeWater_indices.push_back(i);
					closeWater_coords.push_back(waterCoords[i]);
				}
			}

			// If no waters close, we know its completely in the gas phase
			if (closeWater_coords.size() == 0) {
				pinnedClusters.push_back(-1);
				continue;
			}

			// Check other atoms in residue for coordination to test if truly (de)solvated
			bool solvated = false;
			int cluster;
			int numWaters = -1;
			for (const auto& resAtom : titAtom->parent.lock()->atoms) {
				std::vector<float> distances = norm(closeWater_coords, resAtom.lock()->coord);
				int min_dist_index = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));

				// If donor (ie, protonated), use distance to water oxygens
				int water_index;
				if (donor) {
					water_index = closeWater_indices[min_dist_index];
				}
				// Else if deprotonated, use distance to water hydrogens
				else {
					water_index = closeWater_indices[min_dist_index] / 2; // Have to account for 2 hydrogens per water
				}

				// Pin the cluster if within cutoff and has more cluster waters
				int potCluster = clusterIDs[water_index];
				int clusterWater = clusterWaters.at(potCluster);
				if (distances[min_dist_index] < 0.40 && clusterWater > numWaters) {
					cluster = potCluster;
					numWaters = clusterWater;
					solvated = true;
				}
			}

			// If no atoms within 0.35 nm cutoff, assume in gas phase with cluserID -1
			if (!solvated) {
				cluster = -1;
			}

			pinnedClusters.push_back(cluster);
		}

		return pinnedClusters;
	}

	// Finds number of waters in cluster given an inputted list of clusterID's
	std::vector<int> pinClusterWaters(const std::vector<int>& clusterIDs, const std::unordered_map<int, Cluster>& clusters) {
		std::vector<int> clusterWaters;
		for (const int& cluster : clusterIDs) {
			if (clusters.find(cluster) == clusters.end()) {
				throw std::invalid_argument("Cluster ID not found in cluster map");
			}
			else {
				clusterWaters.push_back(clusters.at(cluster).numWaters);
			}
		}

		return clusterWaters;
	}

	// Finds number of a particular ion in a particular cluster
	static std::unordered_map<int, int> pinClusterIons(const std::string& resType, const CoordInfo& coordInfo,
		const std::vector<int>& clusterIDs, const std::unordered_map<int, int>& clusterWaters) {
		
		// Initialize
		std::unordered_map<int, int> clusterIons;
		for (const int& cluster : clusterIDs) {
			clusterIons[cluster] = 0;
		}

		// Get an atom pointer for each ion
		std::vector<std::shared_ptr<Atom>> titAtoms;
		for (const auto& residue : coordInfo.residueMap.at(resType)) {
			titAtoms.push_back(residue->atoms[0].lock());
		}

		// Cluster keeping in mind differences in titratable atoms for H3O+ vs OH-
		std::vector<int> ionClusters;
		if (resType == "HHO") {
			ionClusters = pinCluster(titAtoms, coordInfo.waterOCoords, clusterIDs, clusterWaters, true);
		}
		else if (resType == "OHX") {
			ionClusters = pinCluster(titAtoms, coordInfo.waterHCoords, clusterIDs, clusterWaters, false);
		}
		else {
			std::cerr << "ERROR: Invalid residue type for pinning of cluster ions." << std::endl;
			std::exit(1);
		}

		for (const int& cluster : ionClusters) {
			clusterIons[cluster]++;
		}

		return clusterIons;
	}

	// Calculates pH based on net number of H3O - OH, to ratio of waters
	static float PHcalc(const int& H3O, const int& OH, const int& waters) {
		int netIons = H3O - OH;

		// If net ions positive, use net number of H3O+ 
		if (netIons > 0) {
			return -std::log10(55.0679f * (float)netIons / (float)waters);
		}
		// If net ions negative, use net number of OH- 
		else if (netIons < 0) {
			return 14.0f + std::log10(55.0679f * std::abs((float)OH) / (float)waters);
		}
		// If net ions zero, pH is 7
		else {
			return 7.0f;
		}
	}

	/* Find number of waters, ions, and pH of each cluster(via the ratio of net H3O + &OH - to cluster waters).
	Outputs completed Cluster objects. */
	std::unordered_map<int, Cluster> binClusters(const std::vector<int>& clusterIDs, const CoordInfo& coordInfo, const Config& config) {
		std::unordered_map<int, Cluster> clusterInfo;

		// Water in each cluster
		std::unordered_map<int, int> clusterWaters;
		for (const int& cluster : clusterIDs) {
			clusterWaters[cluster]++;
		}

		// H3O+ and OH- in each cluster
		std::unordered_map<int, int> clusterH3O = pinClusterIons("HHO", coordInfo, clusterIDs, clusterWaters);
		std::unordered_map<int, int> clusterOH = pinClusterIons("OHX", coordInfo, clusterIDs, clusterWaters);

		// pH of each cluster and bin into Cluster objects
		std::unordered_map<int, float> clusterPH;
		for (const int& cluster : clusterIDs) {
			clusterPH[cluster] = PHcalc(clusterH3O[cluster], clusterOH[cluster], clusterWaters[cluster]);
			Cluster clusterObj = { cluster, clusterWaters[cluster], clusterH3O[cluster], clusterOH[cluster], clusterPH[cluster] };
			clusterInfo[cluster] = clusterObj;
		}

		// Cluster -1 represents gas phase 
		clusterInfo[-1] = { -1, 0, 0, 0, 7.0f };

		return clusterInfo;
	}
}