#include "Atmosphere.h"
#include "Atom.h"
#include "Constants.h"
#include "CoordinateManip.h"
#include "Parameters.h"
#include "Vec3D.h"
#include "RNG.h"

#include <cmath>
#include <set>
#include <string>
#include <vector>
#include <random>

int getNumGas(const double boxD) {
	double boxVolume = pow(boxD, 3) * 1e-27; // m^3
	double numDensity = (Constants::standardAtm * Constants::NA) / (Constants::R * Parameters::Get().getGasTemp());
	return static_cast<int>(round(numDensity * boxVolume));
}

unordered_map<string, int> getGasNums() {
	Parameters& params = Parameters::Get();
	static bool humid = params.isHumid();

	int numGas = getNumGas(Parameters::Get().getBoxSize());
	unordered_map<string, int> gases = {
		{"NNN", static_cast<int>(round(numGas * params.getN2Molp() / 100))},
		{"OOO", static_cast<int>(round(numGas * params.getO2Molp() / 100))}
	};

	if (humid)
		gases["SOL"] = static_cast<int>(round(numGas * params.getWaterMolp() / 100));

	return gases;
}

void seedAtmosphere(State& state) {
	Parameters& params = Parameters::Get();
	double boxD = params.getBoxSize();
	int numGas = getNumGas(boxD);

	// Determine grid size
	size_t numPoints = static_cast<size_t>(ceil(cbrt(static_cast<double>(numGas)))) + 1; // Want larger 
	double buffer = 1.0; // Note this must be larger than maxMove
	double step = (boxD - 2 * buffer) / (numPoints - 1);

	// Generate grid
	vector<Vec3D> grid;
	grid.reserve(numPoints * numPoints * numPoints);

	for (size_t ix = 0; ix < numPoints; ++ix) {
		double x = buffer + ix * step;
		for (size_t iy = 0; iy < numPoints; ++iy) {
			double y = buffer + iy * step;
			for (size_t iz = 0; iz < numPoints; ++iz) {
				double z = buffer + iz * step;
				grid.emplace_back(x, y, z);
			}
		}
	}

	// Perturb grid
	double minSpace = 0.6; // (nm), closest distance gases can be, ~2 N2 widths
	double maxMove = (step - minSpace) / 2;
	uniform_real_distribution<double> pd(-maxMove, maxMove);
	for (auto& c : grid) 
		c = c + Vec3D(pd(RNG::gen), pd(RNG::gen), pd(RNG::gen));

	// Remove coordinates that overlap with droplet
	constexpr double inf = numeric_limits<double>::infinity();
	Vec3D maxc(-inf, -inf, -inf);
	Vec3D minc(inf, inf, inf);
	for (const auto& [resType, residues] : state.residueSet) {
		for (const auto& residue : residues) {
			for (const auto& atom : residue->atoms) {
				const Vec3D& coord = atom->coord;
				maxc.x = max(maxc.x, coord.x);
				maxc.y = max(maxc.y, coord.y);
				maxc.z = max(maxc.z, coord.z);
				minc.x = min(minc.x, coord.x);
				minc.y = min(minc.y, coord.y);
				minc.z = min(minc.z, coord.z);
			}
		}
	}

	for (size_t i = grid.size(); i-- > 0;) {
		const Vec3D& c = grid[i];
		if (c.x >= minc.x && c.x <= maxc.x &&
			c.y >= minc.y && c.y <= maxc.y &&
			c.z >= minc.z && c.z <= maxc.z) {
			grid[i] = grid.back();
			grid.pop_back();
		}
	}

	// Ensure grid points equals number of gas molecules
	if (grid.size() > numGas) {
		uniform_int_distribution<int> rd(0, grid.size() - 1);

		set<int, greater<int>> toRemove;
		for (int i = 0; i < grid.size() - numGas; ++i) {
			while (true) {
				int idx = rd(RNG::gen);
				if (!toRemove.contains(idx)) {
					toRemove.insert(idx);
					break;
				}
			}
		}

		for (auto idx : toRemove) {
			grid[idx] = grid.back();
			grid.pop_back();
		}
	}

	// ----- Gas seeding -----
	unordered_map<int, shared_ptr<Residue>> templateResidues = {
		{ 0, parseResidueGRO("n2.gro") }, { 1, parseResidueGRO("o2.gro") }, { 2, parseResidueGRO("sol.gro") }
	};
	vector<double> weights = { params.getN2Molp() / 100, params.getO2Molp() / 100, params.getWaterMolp() / 100 };
	discrete_distribution<int> gd(weights.begin(), weights.end());

	for (size_t i = 0; i < numGas; ++i)
		insertMolec(grid[i], templateResidues.at(gd(RNG::gen)), state);
}