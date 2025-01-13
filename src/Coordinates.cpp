#include "Core.h"
#include <omp.h>
#include <cmath>

namespace Core {

	// Calculate the center of mass of inputted coordinates
	std::array<float, 3> getCOM(const std::vector<std::array<float, 3>>& coords, const int& numCoords) {
		std::array<float, 3> com = { 0.0f, 0.0f, 0.0f };

		for (int i = 0; i <numCoords; i++) {
			for (int j = 0; j < 3; j++) {
				com[j] += coords[i][j];
			}
		}

		for (int j = 0; j < 3; j++) {
			com[j] /= numCoords;
		}

		return com;
	}

	// Calculate the distance between two coordinates
	float getDistance(const std::array<float, 3>& coord1, const std::array<float, 3>& coord2) {
		float s = 0.0f;
		for (int i = 0; i < 3; i++) {
			const float d = coord1[i] - coord2[i];
			s += d * d;
		}
		return std::sqrt(s);
	}

	// Set bond length between atoms by moving atom2 along axis with atom1
	std::array<float, 3> setBondLength(const std::array<float, 3>& coord1, 
		const std::array<float, 3>& coord2, const float& length) {
		std::array<float, 3> direction = subtractVectors(coord2, coord1);
		std::array<float, 3> norm = normalizeVector(direction);
		std::array<float, 3> scaled_direction = scaleVector(norm, length);
		std::array<float, 3> newCoords = addVectors(coord1, scaled_direction);
		return newCoords;
	}

	// Equiavalent of np.linalg.norm, computes distances between point and collection of 3D points
	std::vector<float> norm(const std::vector<std::array<float, 3>>& coords,
		const std::array<float, 3>& refCoord) {

		size_t size = coords.size();
		std::vector<float> distances(size);

		#pragma omp parallel for
		for (size_t i = 0; i < size; ++i) {
			float s = 0.0f;
			for (int k = 0; k < 3; k++) {
				const float d = coords[i][k] - refCoord[k];
				s += d * d;
			}
			distances[i] = std::sqrt(s);
		}

		return distances;
	}

	// Equiavalent of scipy.spatial.distance.cdist, computes distances between two collections of 3D coordinates
	std::vector<std::vector<float>> cdist(const std::vector<std::array<float, 3>>& coords1,
		const std::vector<std::array<float, 3>>& coords2) {

		size_t size1 = coords1.size();
		size_t size2 = coords2.size();
		std::vector<std::vector<float>> distances(size1, std::vector<float>(size2));

		#pragma omp parallel for
		for (size_t i = 0; i < size1; ++i) {
			for (size_t j = 0; j < size2; ++j) {
				float s = 0.0f;
				for (int k = 0; k < 3; k++) {
					const float d = coords1[i][k] - coords2[j][k];
					s += d * d;
				}
				distances[i][j] = std::sqrt(s);
			}
		}

		return distances;
	}
}