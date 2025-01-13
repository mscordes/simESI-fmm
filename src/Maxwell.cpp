#include "Core.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>

namespace Core {

	// Sample n velocities from a Maxwell-Boltzmann distribution in nm/ps
    std::vector<std::array<float, 3>> sampleMaxwell(const int n, const float temperature, const std::string element) {
        float mass {};
        const std::unordered_map<std::string, float> atomMasses = {
            {"H", 1.008f}, {"C", 12.011f}, {"N", 14.007f}, {"O", 15.999f}, {"S", 32.06f}, {"M", 1.000f}
        };
        if (atomMasses.find(element) == atomMasses.end()) {
            std::cerr << "Error: Element " << element << " not found in the database." << std::endl;
			std::exit(1);
        }
        else {
			mass = atomMasses.at(element);
        }
        const float sigma = std::sqrt(0.008314f * temperature / mass); // Standard deviation

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<float> dist(0.0f, sigma);

        // Vector to store the sampled velocities
        std::vector<std::array<float, 3>> velocities;
        velocities.reserve(n);

        // Sample n velocities
        for (int i = 0; i < n; ++i) {
            velocities.emplace_back(std::array<float, 3>{dist(gen), dist(gen), dist(gen)});
        }

        return velocities;
    }
}