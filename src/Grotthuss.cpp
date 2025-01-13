#include "Core.h"
#include <limits>  
#include <algorithm>

namespace Core {

	/* Energy calculations for Grotthuss mechanism requires special consideration, see
		Lars Konermann and Scott Kim, Journal of Chemical Theory and Computation 2022 18 (6), 3781-3794
		DOI: 10.1021/acs.jctc.2c00001 */

	// Potential of single reference point, total Grotthuss energy is delta E = E(donor) - E(acceptor)
	// FIX - this needs to be checked for correspondence
	static float couloumbEnergy(const CoordInfo& coordInfo, const std::vector<float>& charges, const TitratableAtom& refAtom, 
		const TitratableAtom& altAtom, const float& charge) {

		// Get residue objects from atoms
		auto nonWater = refAtom.atom->parent.lock();
		auto water = altAtom.atom->parent.lock();

		// Get coordinates of central atom of each ion (either H3O+ or OH-) as reference point for Couloumb calculation
		std::array<float, 3> refPoint;
		if (refAtom.atom->res_name == "HHO") {
			refPoint =  nonWater->atoms[3].lock()->coord; // Virtual site of H3O+
		}
		else if (refAtom.atom->res_name == "OHX") {
			refPoint = nonWater->atoms[0].lock()->coord; // O of OH-
		}
		else if (refAtom.atom->res_name == "SOL") {
			refPoint = nonWater->atoms[3].lock()->coord; // Virtual site of H2O
		}
		else {
			std::cerr << "Error: Invalid residue type for Grotthuss mechanism." << std::endl;
			exit(1);
		}

		// Avoid using molecules involved in exchange to broadly sample ES environment
		std::vector<int> avoidAtoms;
		for (const auto& atom : nonWater->atoms) {
			avoidAtoms.push_back(atom.lock()->atom_num);
		}
		for (const auto& atom : water->atoms) {
			avoidAtoms.push_back(atom.lock()->atom_num);
		}

		// Begin electrostatics calculation relative to reference point 
		float energy = 0.0f;
		for (int i = 0; i < coordInfo.coordinates.size(); i++) {
			const float dx = coordInfo.coordinates[i][0] - refPoint[0];
			const float dy = coordInfo.coordinates[i][1] - refPoint[1];
			const float dz = coordInfo.coordinates[i][2] - refPoint[2];
			const float distance = std::sqrt(dx * dx + dy * dy + dz * dz);

			// Avoid waters within 0.50 nm of reference point to broadly sample ES environment
			if (distance < 0.50 && coordInfo.atoms[i]->res_name == "SOL") {
				; // Pass
			}
			else if (std::find(avoidAtoms.begin(), avoidAtoms.end(), coordInfo.atoms[i]->atom_num) != avoidAtoms.end()) {
				; // Pass
			}
			// Coulomb contribution
			else {
				energy += (charges[i] * charge) / distance;
			}
		}

		return energy * 138.932f; // Convert to kJ/mol where ke = 138.932 kJ*nm/mol*e^2
	}

	// Compute delta Couloumb energy for Grotthuss mechanism 
	float grotthussEnergy(const CoordInfo& coordInfo, const std::vector<float>& charges, const TitratableAtom& donor, 
		const TitratableAtom& acceptor, const int& hop) {

		float deltaE;

		// Compute Coulomb and apply corrections to ensure proper diffusion coefficient
		if (donor.atom->res_name == "HHO" && acceptor.atom->res_name == "SOL") {
			float charge = 1.0f;
			deltaE = couloumbEnergy(coordInfo, charges, donor, acceptor, charge) - couloumbEnergy(coordInfo, charges, acceptor, donor, charge);
			deltaE -= 25.0f;
		}
		else if (acceptor.atom->res_name == "OHX" && donor.atom->res_name == "SOL") {
			float charge = -1.0f;
			deltaE = couloumbEnergy(coordInfo, charges, donor, acceptor, charge) - couloumbEnergy(coordInfo, charges, acceptor, donor, charge);

			// OH- hops less frequently than H3O+ so apply correction to ensure proper diffusion coefficient
			if (hop > 2) {
				deltaE = 9999.0f;
			}
		}
		else {
			std::cerr << "Error: Invalid residue types for Grotthuss mechanism." << std::endl;
			std::cerr << "Donor: " << donor.atom->res_name << " Acceptor: " << acceptor.atom->res_name << std::endl;
			exit(1);
		}		

		// Enable random comparisons of exchanges with identical energy
		float random = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 0.0001f;

		// Avoid MCMC sampling that used for other exchanges
		if (deltaE < 0.0f) {
			return -9999.f + random;
		}
		else {
			return 9999.f + random;
		}
	}
	
}