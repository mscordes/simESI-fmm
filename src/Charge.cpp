#include "Core.h"
#include <sstream>
#include <unordered_map>
#include <numeric>  
#include <cmath>
#include <string>
#include <cctype>

namespace Core {

	// Partials of all non-protein residues
	std::unordered_map<std::string, std::vector<float>> partialDict = {
		{ "SOL", { 0.00000f, 0.55640f, 0.55640f, -1.11280f }},
		{ "HHO", { 0.00000f, 0.41600f, 0.41600f, -0.24800f, 0.41600f }},
		{ "OHX", { -1.32000f, 0.32000f }},
		{ "NNN", { 0.00000f, 0.00000f }},
		{ "OOO", { 0.00000f, 0.00000f }},
		{ "ATX", { -0.37000f, 0.62000f, 0.09000f, 0.09000f, 0.09000f, -0.76000f, -0.76000f }},
		{ "AHX", { -0.30000f, 0.75000f, 0.09000f, 0.09000f, 0.09000f, -0.55000f, -0.60000f, 0.43000f }},
		{ "NXX", { -1.12500f, 0.37500f, 0.37500f, 0.37500f }},
		{ "NXH", { 0.33000f, -0.32000f, 0.33000f, 0.33000f, 0.33000f }}
	};

	// Get partial charge of each atom in each protein
	std::vector<float> getProtCharges(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms) {
		std::vector<float> protCharges;

		for (const auto& pair : proteinAtoms) {
			std::string chain = pair.first;
			std::string itpName = chain + ".itp";
			std::vector<std::string> itpLines = readFile(itpName);

			for (const auto& line : itpLines) {
				std::vector<std::string> words = splitLine(line);

				if (line.find_first_not_of(' ') != std::string::npos) {
					if (words[1] == "bonds") {
						break;
					}
					else if (words[0][0] != ';' && words[0][0] != '[' && words.size() >= 7) {
						float charge = std::stof(words[6]);
						protCharges.push_back(charge);
					}
				}
			}
		}

		return protCharges;
	}

	// Gets partial charges of all atoms
	std::vector<float> getCharges(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms,
		const std::unordered_map<std::string, int>& numResidues, const std::vector<std::string>& topOrder) {
		
		// Start with protein as first in .top
		std::vector<float> charges = getProtCharges(proteinAtoms);
		for (const auto& resType : topOrder) {
			if (numResidues.find(resType) != numResidues.end() && numResidues.at(resType) > 0) {
				for (int i = 0; i < numResidues.at(resType); i++) {
					charges.insert(charges.end(), partialDict[resType].begin(), partialDict[resType].end());
				}
			}
		}

		return charges;
	}

	// Sum a vector of charges
	int getNetCharge(const std::vector<float>& charges) {
		int netCharge = std::round(std::accumulate(charges.begin(), charges.end(), 0.0f));
		return netCharge;
	}
}