#include "CoordinateManip.h"
#include "Atom.h"
#include "Constants.h"
#include "MathUtils.h"
#include "Protein.h"
#include "Vec3D.h"
#include "RNG.h"
#include "StringUtils.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>

vector<Vec3D> getAllCoords(const State& state) {
	vector<Vec3D> coords;
	for (const auto& [_, residues] : state.residueSet)
		for (const auto& r : residues)
			for (const auto& a : r->atoms)
				coords.push_back(a->coord);
	return coords;
}

double getProteinMass(const State& state) {
	double mass = 0;
	for (const auto& [chain, protein] : state.proteinMap)
		mass += protein->getMass();
	return mass;
}

vector<Vec3D> getProteinCoords(const State& state) {
	size_t count = 0;
	for (const auto& [_, p] : state.proteinMap)
		for (const auto& r : p->residues)
			count += r->atoms.size();

	vector<Vec3D> coords;
	coords.reserve(count);
	for (const auto& [chain, p] : state.proteinMap)
		for (const auto& r : p->residues) 
			for (const auto& a : r->atoms) 
				coords.push_back(a->coord);
	return coords;
}

vector<Vec3D> getProteinCarbonCoords(const State& state) {
	size_t count = 0;
	for (const auto& [_, p] : state.proteinMap)
		for (const auto& r : p->residues)
			for (const auto& a : r->atoms)
				if (a->element == 'C')
					count += 1;

	vector<Vec3D> coords;
	coords.reserve(count);
	for (const auto& [chain, p] : state.proteinMap)
		for (const auto& r : p->residues)
			for (const auto& a : r->atoms)
				if (a->element == 'C')
					coords.push_back(a->coord);
	return coords;
}

vector<Vec3D> getProteinSkeletonCoords(const State& state) {
	size_t count = 0;
	for (const auto& [_, p] : state.proteinMap)
		for (const auto& r : p->residues)
			count += 1;

	vector<Vec3D> coords;
	coords.reserve(count);
	for (const auto& [chain, p] : state.proteinMap)
		for (const auto& r : p->residues)
			coords.push_back(r->atoms[0]->coord);
	return coords;
}


shared_ptr<Atom> findClosestWaterO(double& dist, const Vec3D& ref, const WaterInfo& waterInfo) { 
	double min_d2 = numeric_limits<double>::infinity();
	size_t min_idx = numeric_limits<size_t>::max();
	const double rx = ref.x, ry = ref.y, rz = ref.z;
	const auto& coords = waterInfo.Ocoords;

	for (size_t i = 0, n = coords.size(); i < n; ++i) {
		const Vec3D& c = coords[i];
		const double dx = rx - c.x;
		const double dy = ry - c.y;
		const double dz = rz - c.z;
		const double d2 = dx * dx + dy * dy + dz * dz;
		if (d2 < min_d2) {
			min_d2 = d2;
			min_idx = i;
		}
	}

	if (min_idx == numeric_limits<size_t>::max())
		throw runtime_error("No waters in State when finding closest water O.");

	dist = sqrt(min_d2);
	return waterInfo.waters[min_idx]->atoms[0];
}

shared_ptr<Atom> findClosestWaterH(double& dist, const Vec3D& ref, const WaterInfo& waterInfo) {
	double min_d2 = numeric_limits<double>::infinity();
	size_t min_idx = numeric_limits<size_t>::max();
	const double rx = ref.x, ry = ref.y, rz = ref.z;
	const auto& coords = waterInfo.Hcoords;

	for (size_t i = 0, n = coords.size(); i < n; ++i) {
		const Vec3D& c = coords[i];
		const double dx = rx - c.x;
		const double dy = ry - c.y;
		const double dz = rz - c.z;
		const double d2 = dx * dx + dy * dy + dz * dz;
		if (d2 < min_d2) {
			min_d2 = d2;
			min_idx = i;
		}
	}

	if (min_idx == numeric_limits<size_t>::max())
		throw runtime_error("No waters in State when finding closest water H.");

	dist = sqrt(min_d2);
	return waterInfo.waters[min_idx / 2]->atoms[(min_idx % 2) + 1]; // 2 H's per water
}

void centerState(State& state, const Vec3D& boxD) {
	Vec3D center = boxD / 2;
	const auto pCoords = getProteinCoords(state);
	Vec3D pCenter = getCenter(pCoords);

	Vec3D translation = center - pCenter;
	for (const auto& [resType, residues] : state.residueSet) 
		for (const auto& r : residues) 
			for (const auto& a : r->atoms) 
				a->coord = a->coord + translation;
}

void centerState(State& state, const double boxScalar) {
	centerState(state, { boxScalar, boxScalar, boxScalar });
}

// Make as small a cubic box as possible that can solvate a given protein + droplet
Vec3D centerSmallBox(State& state) {
	constexpr double inf = numeric_limits<double>::infinity();
	Vec3D maxc(-inf, -inf, -inf);
	Vec3D minc(inf, inf, inf);
	for (const auto& [_, p] : state.proteinMap) {
		for (const auto& r : p->residues) {
			for (const auto& a : r->atoms) {
				const Vec3D& coord = a->coord;
				maxc.x = max(maxc.x, coord.x);
				maxc.y = max(maxc.y, coord.y);
				maxc.z = max(maxc.z, coord.z);
				minc.x = min(minc.x, coord.x);
				minc.y = min(minc.y, coord.y);
				minc.z = min(minc.z, coord.z);
			}
		}
	}

	Vec3D boxD = (maxc - minc);
	boxD = boxD + (Parameters::Get().getDropletSize() * 2); // Larger than protein + droplet;
	centerState(state, boxD);
	return boxD;
}

void carveDroplet(State& state, const double envelope) {
	auto it = state.residueSet.find("SOL");
	if (it == state.residueSet.end()) return;
	auto& waters = it->second;

	vector<Vec3D> skeleton;
	vector<Vec3D> carbons;
	for (const auto& [_, p] : state.proteinMap) {
		for (const auto& r : p->residues) {
			skeleton.push_back(r->atoms[0]->coord);
			for (const auto& a : r->atoms)
				if (a->element == 'C')
					carbons.push_back(a->coord);
		}
	}

	const double env_2 = envelope * envelope;
	const double skel_2 = 2 * envelope * envelope; 

	unordered_set<shared_ptr<Residue>> toDelete;
	toDelete.reserve(waters.size());

	for (const auto& w : waters) {
		const Vec3D& O = w->atoms[0]->coord;
		bool inside = false;

		for (const Vec3D& c : skeleton) {
			if (O.distance_sq(c) < skel_2) { inside = true; break; }
		}

		if (!inside) {
			for (const Vec3D& c : carbons) {
				if (O.distance_sq(c) < env_2) { inside = true; break; }
			}
		}

		if (!inside)
			toDelete.insert(w);
	}

	for (const auto& w : toDelete)
		waters.erase(w);
}

// Parse a .gro file of a single residue
shared_ptr<Residue> parseResidueGRO(const string& groFile) {
	ifstream file(groFile);
	if (!file.is_open()) throw runtime_error("Error opening .gro file: " + groFile);

	// Skip the first two lines
	string line;
	getline(file, line);
	getline(file, line);

	shared_ptr<Residue> r;
	bool initialized = false;

	while (getline(file, line)) {

		// Ignore last line with box dimensions
		if (line.size() < 44) break;

		if (!initialized) {
			string name = trim(line.substr(5, 3));
			r = make_shared<Residue>("-1" + name, name, -1, 'Z');

			if (name != "SOL" && name != "NNN" && name != "OOO") { // SOL can be either donor/acceptor so handle specially later
				if (auto p = Constants::defaultPkaVals.find(name); p != Constants::defaultPkaVals.end())
					r->pka = p->second;
				else throw runtime_error("Missing defaultPkaVals for " + name);

				if (auto g = Constants::GPBs.find(name); g != Constants::GPBs.end())
					r->gpb = g->second;
				else throw runtime_error("Missing GPBs for " + name);
			}

			initialized = true;
		}

		string atomName = trim(line.substr(9, 6));
		r->addAtom(
			make_shared<Atom>(
				atomName,
				atomName[0],
				Vec3D(
					stod(line.substr(20, 8)),
					stod(line.substr(28, 8)),
					stod(line.substr(36, 8))),
				Vec3D(0, 0, 0),
				r
			)
		);
	}
	file.close();

	return r;
}

shared_ptr<Residue> cloneResidue(const Residue& src) {
	auto r = make_shared<Residue>(src.ID, src.name, src.num, src.chain);
	r->gpb = src.gpb;
	r->pka = src.pka;

	r->atoms.reserve(src.atoms.size());
	for (const auto& a : src.atoms) {
		auto newAtom = make_shared<Atom>(*a);
		newAtom->parent = r;
		r->atoms.push_back(newAtom);
	}
	return r;
}

// Insert a molecule of a given type into the droplet
void insertMolecIntoDroplet(
	const string& resType,			
	const unordered_map<string, int>& tooSeed,
	State& state,					 
	const vector<Vec3D>& protCoords, 
	vector<Vec3D>& coords,			 
	const Vec3D& boxD,					
	double minSep)				
{
	int n;
	if (auto it = tooSeed.find(resType); it != tooSeed.end()) n = it->second;
	else return;

	auto templateResidue = parseResidueGRO(Constants::templateGroNames.at(resType));
	const double envelope = Parameters::Get().getDropletSize();
	const double envelope_sq = envelope * envelope;
	const double minSep_sq = minSep * minSep;
	const double wallBuffer = 0.5; // nm

	uniform_real_distribution<double> dx(wallBuffer, boxD.x - wallBuffer);
	uniform_real_distribution<double> dy(wallBuffer, boxD.y - wallBuffer);
	uniform_real_distribution<double> dz(wallBuffer, boxD.z - wallBuffer);

	for (int i = 0; i < n; ++i) {
		int attempts = 0;
		bool placed = false;

		while (!placed) {
			if (++attempts > 10'000)
				throw runtime_error("Could not insert molecule into droplet. Try increasing droplet size.");

			Vec3D newLoc(dx(RNG::gen), dy(RNG::gen), dz(RNG::gen));

			// Must be within envelope of protein
			bool insideEnvelope = false;
			for (const auto& pc : protCoords) {
				if (newLoc.distance_sq(pc) < envelope_sq) { insideEnvelope = true; break; }
			}
			if (!insideEnvelope) continue;

			auto newResidue = cloneResidue(*templateResidue);

			// Check if the molecule is too close to another molecule
			bool tooClose = false;
			for (const auto& atom : newResidue->atoms) {
				Vec3D& newCoord = atom->coord;
				newCoord = newCoord + newLoc;

				for (const auto& c : coords) {
					if (newCoord.distance_sq(c) < minSep_sq) {
						tooClose = true;
						break;
					}
				}
				if (tooClose) break;
			}
			if (tooClose) continue;

			// Update coords
			for (const auto& atom : newResidue->atoms)
				coords.push_back(atom->coord);

			state.residueSet[newResidue->name].insert(newResidue);
			placed = true;
		}
	}
}

// Insert molecule at inputted location given a template
void insertMolec(const Vec3D& newLoc, const shared_ptr<Residue>& templateResidue, State& state) {
	auto newResidue = cloneResidue(*templateResidue);
	for (const auto& atom : newResidue->atoms)
		atom->coord = atom->coord + newLoc;
	state.residueSet[newResidue->name].insert(newResidue);
}

// Fix annoying bug where will pdb2gmx will not recognize disulfide bond if bond extends beyond 2.0 ± 0.2 A
// Return true if coordinates had to be modified
void fixDisulfides(State& state, RunFlags& flags) {
	auto it = state.residueSet.find("CYS");
	if (it == state.residueSet.end()) return;

	const auto& cysteines = it->second;
	auto findSG = [](const shared_ptr<Residue>& r) -> shared_ptr<Atom> {
		for (auto& a : r->atoms) if (a->name == "SG") return a;
		return nullptr;
		};

	for (const auto& cys1 : cysteines) {
		auto sg1 = findSG(cys1);
		if (!sg1) continue;

		for (const auto& cys2 : cysteines) {
			if (cys1.get() == cys2.get()) continue;

			auto sg2 = findSG(cys2);
			if (!sg2) continue;

			const auto& c1 = sg1->coord;
			const auto& c2 = sg2->coord;
			double d = c1.distance(c2);
			if (d > 0.25) continue;

			if (d < 0.19 || d > 0.21) {
				flags.createRunFile = true;

				double targetLength = 0.200; 
				Vec3D dir = (c2 - c1).normalize();
				Vec3D midpoint = (c1 + c2) / 2.0;

				sg1->coord = midpoint - dir * targetLength / 2.0;
				sg2->coord = midpoint + dir * targetLength / 2.0;
			}
		}
	}
}

WaterInfo getWaterInfo(const State& state) {
	auto it = state.residueSet.find("SOL");
	if (it == state.residueSet.end()) return {};
	const auto& waters = it->second;

	WaterInfo info;
	info.waters.reserve(waters.size());
	info.Ocoords.reserve(waters.size());
	info.Hcoords.reserve(waters.size() * 2);
	info.coords.reserve(waters.size());

	for (const auto& water : waters) {    
		const auto& atoms = water->atoms;
		info.waters.push_back(water);
		info.Ocoords.push_back(atoms[0]->coord);
		info.Hcoords.push_back(atoms[1]->coord);
		info.Hcoords.push_back(atoms[2]->coord);
		info.coords.push_back({ atoms[0]->coord,
								atoms[1]->coord,
								atoms[2]->coord, 
								atoms[3]->coord});
	}
	return info;
}

// Gets partial charges of all atoms excluding protein, water, and atmosphere
unordered_map<string, vector<Vec3D>> getSoluteCoords(const State& state) {
	unordered_map<string, vector<Vec3D>> coords;
	coords.reserve(Constants::topOrder.size()); 

	for (const auto& resType : Constants::topOrder) {
		if (resType == "SOL" || resType == "NNN" || resType == "OOO") continue;

		auto it = state.residueSet.find(resType);
		if (it == state.residueSet.end()) continue;

		const auto& residues = it->second;
		if (residues.empty()) continue;

		vector<Vec3D> sCoords;
		size_t numResAtoms = (*residues.begin())->atoms.size();
		sCoords.reserve(residues.size() * numResAtoms);

		for (const auto& r : residues)
			for (const auto& a : r->atoms)
				sCoords.push_back(a->coord);

		coords.emplace(resType, move(sCoords));
	}
	return coords;
}

unordered_set<shared_ptr<Residue>> getGaseousWaters(const State& state) {
	auto it = state.residueSet.find("SOL");
	if (it == state.residueSet.end()) return {};

	const auto& waters = it->second;
	if (waters.size() == 0) return {};

	const auto& skeleton = getProteinSkeletonCoords(state);
	const double cutoff_sq = 5.0 * 5.0; // 5 nm cutoff for now

	unordered_set<shared_ptr<Residue>> gaseousWaters;
	for (const auto& w : waters) {
		const auto& wc = w->atoms[0]->coord;

		bool isGas = true;
		for (const auto& c : skeleton) {
			if (wc.distance_sq(c) < cutoff_sq) {
				isGas = false;
				break;
			}
		}
		if (isGas) gaseousWaters.insert(w);
	}

	return gaseousWaters;
}