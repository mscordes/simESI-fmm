#include "Core.h"
#include <cmath>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>

namespace Core {

	// Constants
	const float kB = 8314.0f; // Boltzmann constant (J/K)

	// Helper function for the error function approximation
	static float erf(float x) {
		// Use a numerical approximation for the error function
		// Abramowitz and Stegun formula
		const float a1 = 0.254829592;
		const float a2 = -0.284496736;
		const float a3 = 1.421413741;
		const float a4 = -1.453152027;
		const float a5 = 1.061405429;
		const float p = 0.3275911;

		int sign = (x < 0) ? -1 : 1;
		x = std::fabs(x);

		float t = 1.0 / (1.0 + p * x);
		float y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

		return sign * y;
	}

	// Maxwell-Boltzmann CDF function
	static float maxwell_CDF(float v, float T, float m) {
		// Precompute constants
		float sqrt_term = std::sqrt(m / (2 * kB * T));
		float prefactor = std::sqrt(2.0 * m / (M_PI * kB * T));

		// Compute the CDF
		float erf_part = erf(sqrt_term * v);
		float exp_part = prefactor * v * std::exp(-m * v * v / (2 * kB * T));

		return erf_part - exp_part;
	}

	// Generate the Maxwell-Boltzmann CDF over a range of speeds
	static std::vector<std::pair<float, float>> compute_CDF(float T, float m, float max_speed, int num_points) {
		std::vector<std::pair<float, float>> cdf_values;
		float step = max_speed / num_points;

		for (int i = 0; i <= num_points; ++i) {
			float v = i * step;
			float cdf = maxwell_CDF(v, T, m);
			cdf_values.emplace_back(v, cdf);
		}

		return cdf_values;
	}

	// Returns the ECDF of a given sample
	static std::vector<std::pair<float, float>> compute_ECDF(const std::vector<float>& data) {
		std::vector<float> sorted_data = data;
		std::sort(sorted_data.begin(), sorted_data.end());

		std::vector<std::pair<float, float>> ecdf; // Pair of (value, ECDF)
		size_t n = sorted_data.size();

		for (size_t i = 0; i < n; ++i) {
			float value = sorted_data[i];
			float fraction = static_cast<float>(i + 1) / n; // ECDF value
			ecdf.emplace_back(value, fraction);
		}

		return ecdf;
	}

	// Function to calculate the KS statistic
	static float ksStatistic(const std::vector<std::pair<float, float>>& empirical_ecdf,
		const std::vector<std::pair<float, float>>& actual_ecdf) {
		float max_diff = 0.0f;

		// Pointers to traverse both ECDFs
		size_t i = 0, j = 0;
		size_t n_empirical = empirical_ecdf.size();
		size_t n_actual = actual_ecdf.size();

		while (i < n_empirical && j < n_actual) {
			// Compare ECDF values at the closest points
			float x_empirical = empirical_ecdf[i].first;
			float x_actual = actual_ecdf[j].first;

			float value_empirical = empirical_ecdf[i].second;
			float value_actual = actual_ecdf[j].second;

			// Compute the absolute difference between ECDFs
			max_diff = std::max(max_diff, std::fabs(value_empirical - value_actual));

			// Increment the pointer for the smaller x-value
			if (x_empirical < x_actual) {
				++i;
			}
			else if (x_actual < x_empirical) {
				++j;
			}
			else {
				++i;
				++j;
			}
		}

		return max_diff;
	}

	// Gas phase velocities can become distorted from stitching. Compare to Maxwell boltzmann distribution, 
	// if deviated, correct gas vels. Modifies coordInfo.atoms in place.
	bool correctGasVelocities(std::vector<std::shared_ptr<Atom>>& atoms, float temperature) {
		bool corrected = false;

		// Each gas has its own Maxwell-Boltzmann distribution
		std::vector<std::string> gasTypes = { "NNN", "OOO" };
		for (const auto& gasType : gasTypes) {

			// Get actual speeds
			std::vector<float> actualSpeeds;
			for (std::vector<std::shared_ptr<Atom>>::reverse_iterator riter = atoms.rbegin(); riter != atoms.rend(); ++riter) {
				if (std::find(gasTypes.begin(), gasTypes.end(), (*riter)->res_name) != gasTypes.end()) {
					if ((*riter)->res_name == gasType) {
						float speed = std::hypot((*riter)->velocity[0], (*riter)->velocity[1], (*riter)->velocity[2]);
						actualSpeeds.push_back(speed * 1000); // Convert to m/s
					}
				}
				else {
					break;
				}
			}

			// Compute ECDF of the actual speeds
			std::vector<std::pair<float, float>> empirical_ecdf = compute_ECDF(actualSpeeds);

			// Compute the  CDF 
			float max_speed = *std::max_element(actualSpeeds.begin(), actualSpeeds.end());
			float mass = gasType == "NNN" ? 14.007f : 15.999f;
			std::vector<std::pair<float, float>> actual_cdf = compute_CDF(temperature, mass, max_speed, actualSpeeds.size());

			// Compute the KS statistic 
			float ks_stat = ksStatistic(empirical_ecdf, actual_cdf);

			// Actual gas vels don't conform to Maxwell-Boltzmann distribution given bonded forces
			// Use KS stat instead, we're really looking for massive deviations due to stitching
			if (ks_stat > 0.2) {
				corrected = true;

				// Generate new velocities from the Maxwell-Boltzmann distribution for N2 and O2
				std::string element = gasType.substr(0, 1);
				std::vector<std::array<float, 3>> expectedVelocities = sampleMaxwell(actualSpeeds.size(), temperature, element);
				
				// Set the new O2 and N2 velocities
				int count = 0;
				for (std::vector<std::shared_ptr<Atom>>::reverse_iterator riter = atoms.rbegin(); riter != atoms.rend(); ++riter) {
					if (std::find(gasTypes.begin(), gasTypes.end(), (*riter)->res_name) != gasTypes.end()) {
						if ((*riter)->res_name == gasType) {
							(*riter)->velocity = expectedVelocities[count];
							count++;
						}
					}
					else {
						break;
					}
				}
				if (count != actualSpeeds.size()) {
					std::cerr << "Error: Incorrect number of gas atoms found when regenerating velocities." << std::endl;
					std::exit(1);
				}
			}
		}

		return corrected;
	}
	
	// Ammonia molecules don't play nice in the gas phase so force them to behave
	void fixAmmoniaGas(const Config& config, CoordInfo& coordInfo, bool& gasCorrect) {
		if (coordInfo.residueMap["NXX"].empty() || config.nh3_vapor < 0.0001) {
			return;
		}

		std::vector<std::shared_ptr<Atom>> refAtoms = readGROatoms("nh3.gro");
		std::vector<std::array<float, 3>> refCoords = extractCoordinates(refAtoms);

		// Parse angle of ammonia H-H-H to see if distorted, correct if bad
		for (auto& residue : coordInfo.residueMap["NXX"]) {

			// Check that residue is gas phase, if not skip
			bool inGas = true;
			std::array<float, 3> loc = residue->atoms[0].lock()->coord;
			for (const auto& coord : coordInfo.waterOCoords) {
				const float dx = loc[0] - coord[0];
				const float dy = loc[1] - coord[1];
				const float dz = loc[2] - coord[2];
				const float distanceSquared = dx * dx + dy * dy + dz * dz;

				if (distanceSquared < 0.25) {
					inGas = false;
					break;
				}
			}
			if (inGas) {
				for (const auto& coord : coordInfo.proteinCoords) {
					const float dx = loc[0] - coord[0];
					const float dy = loc[1] - coord[1];
					const float dz = loc[2] - coord[2];
					const float distanceSquared = dx * dx + dy * dy + dz * dz;

					if (distanceSquared < 0.25) {
						inGas = false;
						break;
					}
				}
			}
			if (!inGas) {
				continue;
			}

			// Check angle and reset molecule if necessary
			float angle = getAngle(residue->atoms[1].lock()->coord, residue->atoms[2].lock()->coord, residue->atoms[3].lock()->coord);
			if (angle > 125.0f || angle < 93.0) {
				gasCorrect = true;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 3; j++) {
						residue->atoms[i].lock()->coord[j] = refCoords[i][j] + loc[j];
					}
					std::string element = residue->atoms[i].lock()->element;
					if (element == "H") {
						residue->atoms[i].lock()->velocity = sampleMaxwell(1, config.gas_temp, element)[0];
					}
				}
			}
		}
	}

	// Acetic acid molecules don't play nice in the gas phase so force them to behave
	void fixAceticGas(const Config& config, CoordInfo& coordInfo, bool& gasCorrect) {
		if (coordInfo.residueMap["AHX"].empty() || config.ach_vapor < 0.0001) {
			return;
		}

		std::vector<std::shared_ptr<Atom>> refAtoms = readGROatoms("aceh.gro");
		std::vector<std::array<float, 3>> refCoords = extractCoordinates(refAtoms);

		// Parse velocities as these become high as ach begins to shake itself apart
		for (auto& residue : coordInfo.residueMap["AHX"]) {

			// Check that residue is gas phase, if not skip
			bool inGas = true;
			std::array<float, 3> loc = residue->atoms[0].lock()->coord;
			for (const auto& coord : coordInfo.waterOCoords) {
				const float dx = loc[0] - coord[0];
				const float dy = loc[1] - coord[1];
				const float dz = loc[2] - coord[2];
				const float distanceSquared = dx * dx + dy * dy + dz * dz;

				if (distanceSquared < 0.25) {
					inGas = false;
					break;
				}
			}
			if (inGas) {
				for (const auto& coord : coordInfo.proteinCoords) {
					const float dx = loc[0] - coord[0];
					const float dy = loc[1] - coord[1];
					const float dz = loc[2] - coord[2];
					const float distanceSquared = dx * dx + dy * dy + dz * dz;

					if (distanceSquared < 0.25) {
						inGas = false;
						break;
					}
				}
			}
			if (!inGas) {
				continue;
			}

			// Check vels
			bool correct = false;
			bool exit_loops = false;
			for (const auto& weak_atom : residue->atoms) {
				auto atom = weak_atom.lock();
				for (int i = 0; i < 3; i++) {
					if (std::abs(atom->velocity[i]) > 7.0f) {
						correct = true;
						exit_loops = true;
						break;
					}
				}
				if (exit_loops) break;
			}

			// Reset molecule if necessary
			if (correct) {
				gasCorrect = true;
				for (int i = 0; i < refCoords.size(); i++) {
					for (int j = 0; j < 3; j++) {
						residue->atoms[i].lock()->coord[j] = refCoords[i][j] + loc[j];
					}
					std::string element = residue->atoms[i].lock()->element;
					residue->atoms[i].lock()->velocity = sampleMaxwell(1, config.gas_temp, element)[0];
				}
			}
		}
	}

}