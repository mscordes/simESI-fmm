#include "DoExchanges.h"
#include "Constants.h"
#include "MathUtils.h"
#include "Transforms.h"

#include <utility>
 
static void transferCoordVel(shared_ptr<Atom>& a1, shared_ptr<Atom>& a2) {
	a2->coord = a1->coord;
	a2->velocity = a1->velocity;
}

static void swapCoordVel(shared_ptr<Atom>& a1, shared_ptr<Atom>& a2) {
	swap(a1->coord, a2->coord);
	swap(a1->velocity, a2->velocity);
}

static void swapAtoms(shared_ptr<Atom>& a1, shared_ptr<Atom>& a2) {
	swap(a1->name, a2->name);
	swap(a1->element, a2->element);
	swap(a1->coord, a2->coord);
	swap(a1->velocity, a2->velocity);
}

// Grotthuss mechanism for H3O+, given that products = reactants can just swap the two
static void doGrotthussH3O(
	const shared_ptr<Atom>& h, 
	const shared_ptr<Atom>& a, 
	shared_ptr<Residue>& rh,
	shared_ptr<Residue>& ra,
	State& state)
{
	for (int i = 0; i < 4; ++i)
		swapCoordVel(rh->atoms[i], rh->atoms[i]);

	if (h->name == "HW1")		
		swapCoordVel(rh->atoms[1], rh->atoms[4]);
	else if (h->name == "HW2")	
		swapCoordVel(rh->atoms[2], rh->atoms[4]);

	hydroniumTransform(rh);
	waterTransform(ra);
}

// Grotthuss mechanism for OH-, given that products = reactants can just swap the two
static void doGrotthussOH(
	const shared_ptr<Atom>& h,
	const shared_ptr<Atom>& a,
	shared_ptr<Residue>& rh,
	shared_ptr<Residue>& ra,
	State& state)
{
	swapCoordVel(rh->atoms[0], ra->atoms[0]);
	if (h->name == "HW2")
		swapCoordVel(rh->atoms[2], ra->atoms[1]);
	else if (h->name == "HW3") 
		swapCoordVel(rh->atoms[1], ra->atoms[1]);

	waterTransform(rh);
}

static void createHydroxide(
	const shared_ptr<Atom>& h,
	shared_ptr<Residue>& rh,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "OHX";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	// O
	transferCoordVel(rh->atoms[0], product->atoms[0]); 

	// H
	if (h->name == "HW2") 
		transferCoordVel(rh->atoms[2], product->atoms[1]);
	else if (h->name == "HW3")
		transferCoordVel(rh->atoms[1], product->atoms[1]);
	else
		throw runtime_error("Unexpected hydrogen name during createHydroxide(): " + h->name);

	int clusterID = rh->cluster;
	product->cluster = clusterID;
	auto& cluster = clusters.at(clusterID);
	if (!cluster.water || !cluster.oh)
		throw runtime_error("createHydroxide: Cluster water, or OH- uninitialized.");
	*cluster.water -= 1;
	*cluster.oh += 1;
	cluster.calculatePH();

	toRemove[rh->name].push_back(rh);
	state.residueSet[string(resType)].insert(product);
}

static void createWaterFromHydronium(
	const shared_ptr<Atom>& h,
	shared_ptr<Residue>& rh,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "SOL";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	// O and M
	transferCoordVel(rh->atoms[0], product->atoms[0]);
	transferCoordVel(rh->atoms[3], product->atoms[3]);

	// H's
	if (h->name == "HW1") { 
		transferCoordVel(rh->atoms[2], product->atoms[1]);
		transferCoordVel(rh->atoms[4], product->atoms[2]);
	}
	else if (h->name == "HW2") {
		transferCoordVel(rh->atoms[1], product->atoms[1]);
		transferCoordVel(rh->atoms[4], product->atoms[2]);
	}
	else if (h->name == "HW3") {
		transferCoordVel(rh->atoms[1], product->atoms[1]);
		transferCoordVel(rh->atoms[2], product->atoms[2]);
	}
	else
		throw runtime_error("Unexpected hydrogen name during createWaterFromHydronium(): " + h->name);

	waterTransform(product);

	int clusterID = rh->cluster;
	product->cluster = clusterID;
	auto& cluster = clusters.at(clusterID);
	if (!cluster.water || !cluster.h3o)
		throw runtime_error("createWaterFromHydronium: Cluster water, or H3O+ uninitialized.");
	*cluster.water += 1;
	*cluster.h3o -= 1;
	cluster.calculatePH();

	toRemove[rh->name].push_back(rh);
	state.residueSet[string(resType)].insert(product);
}

static void createAcetate(
	const shared_ptr<Atom>& h,
	shared_ptr<Residue>& rh,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "ATX";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	for (size_t i = 0; i < product->atoms.size(); ++i)
		transferCoordVel(rh->atoms[i], product->atoms[i]);

	product->cluster = rh->cluster;

	toRemove[rh->name].push_back(rh);
	state.residueSet[string(resType)].insert(product);
}

static void createAmmonia(
	const shared_ptr<Atom>& h,
	shared_ptr<Residue>& rh,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "NXX";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	auto& rh_atoms = rh->atoms;
	auto& p_atoms = product->atoms;
	transferCoordVel(rh_atoms[1], p_atoms[0]);

	if (h->name == "HZ1") {
		transferCoordVel(rh_atoms[2], p_atoms[1]);
		transferCoordVel(rh_atoms[3], p_atoms[2]);
		transferCoordVel(rh_atoms[4], p_atoms[3]);
	}
	else if (h->name == "HZ2") {
		transferCoordVel(rh_atoms[0], p_atoms[1]);
		transferCoordVel(rh_atoms[3], p_atoms[2]);
		transferCoordVel(rh_atoms[4], p_atoms[3]);
	}
	else if (h->name == "HZ3") {
		transferCoordVel(rh_atoms[0], p_atoms[1]);
		transferCoordVel(rh_atoms[2], p_atoms[2]);
		transferCoordVel(rh_atoms[4], p_atoms[3]);
	}
	else if (h->name == "HZ4") {
		transferCoordVel(rh_atoms[0], p_atoms[1]);
		transferCoordVel(rh_atoms[2], p_atoms[2]);
		transferCoordVel(rh_atoms[3], p_atoms[3]);
	}
	else
		throw runtime_error("Bad donor H name of " + h->name + " for NXH.");

	product->cluster = rh->cluster;

	toRemove[rh->name].push_back(rh);
	state.residueSet[string(resType)].insert(product);
}

static void deprotonateAminoAcid(
	const shared_ptr<Exchange>& e,
	const shared_ptr<Atom>& h,
	shared_ptr<Residue>& rh)
{
	auto& atoms = rh->atoms;
	string atomToErase;

	// Correct Residue atom ordering/naming, update acceptors, and mark atom for deletion
	if (e->isNterm) {
		atomToErase = "H3";
		shared_ptr<Atom> h1, h2, h3;

		for (auto& a : atoms) {
			const auto& aName = a->name;
			if		(aName == "H1") h1 = a;
			else if (aName == "H2") h2 = a;
			else if (aName == "H3") h3 = a;
		}

		if (h->name == "H1") swapCoordVel(h1, h3);
		else if (h->name == "H2") swapCoordVel(h2, h3);
	}
	else if (e->isCterm) {
		atomToErase = "HT2";
	}
	else if (rh->name == "LYS") {
		atomToErase = "HZ3";
		shared_ptr<Atom> hz1, hz2, hz3;

		for (auto& a : atoms) {
			const auto& aName = a->name;
			if		(aName == "HZ1") hz1 = a;
			else if (aName == "HZ2") hz2 = a;
			else if (aName == "HZ3") hz3 = a;
		}

		if (h->name == "HZ1") swapCoordVel(hz1, hz3);
		else if (h->name == "HZ2") swapCoordVel(hz2, hz3);
	}
	else if (rh->name == "ARG") {
		atomToErase = "HH12";
		shared_ptr<Atom> nh1, nh2, hh11, hh12, hh21, hh22;

		for (auto& a : atoms) {
			const auto& aName = a->name;
			if		(aName == "NH1") nh1 = a;
			else if (aName == "NH2") nh2 = a;
			else if (aName == "HH11") hh11 = a;
			else if (aName == "HH12") hh12 = a;
			else if (aName == "HH21") hh21 = a;
			else if (aName == "HH22") hh22 = a;
		}

		if (h->name == "HH11") {
			swapCoordVel(hh11, hh12);
		}
		else if (h->name == "HH21") {
			swapCoordVel(nh1, nh2);
			swapCoordVel(hh11, hh22);
			swapCoordVel(hh12, hh21);
		}
		else if (h->name == "HH22") {
			swapCoordVel(nh1, nh2);
			swapCoordVel(hh11, hh21);
			swapCoordVel(hh12, hh22);
		}
	}
	else if (rh->name == "ASP" || rh->name == "GLU") {
		atomToErase = (rh->name == "ASP") ? "HD2" : "HE2";
	}
	else if (rh->name == "HIS") {
		atomToErase = h->name;
		string accName;
		auto& hatoms = rh->atoms;

		// CHARMM changes ordering of HISD/HISE relative to HISP
		if (h->name == "HD1") {
			atomToErase = h->name;
			accName = "ND1";
			swapAtoms( hatoms[7], hatoms[12]);
			swapAtoms( hatoms[8],  hatoms[9]);
			swapAtoms( hatoms[9], hatoms[14]);
			swapAtoms(hatoms[10], hatoms[15]);
			swapAtoms(hatoms[11], hatoms[15]);
			swapAtoms(hatoms[12], hatoms[15]);
			swapAtoms(hatoms[14], hatoms[15]);
		}
		else if (h->name == "HE2") {
			atomToErase = h->name;
			accName = "NE2";
			swapAtoms( hatoms[7], hatoms[12]);
			swapAtoms( hatoms[8], hatoms[13]);
			swapAtoms(hatoms[10], hatoms[14]);
			swapAtoms(hatoms[12], hatoms[15]);
			swapAtoms(hatoms[13], hatoms[14]);
			swapAtoms(hatoms[14], hatoms[15]);
		}
		else
			throw runtime_error("Improper atom type of " + h->name + " for HIS deprotonation.");
	}
	else 
		throw runtime_error("Could not identify protein residue type " + rh->name + " for deprotonation.");


	// Erase donor hydrogen
	size_t idx = atoms.size();
	for (size_t i = 0; i < atoms.size(); i++) {
		if (atoms[i]->name == atomToErase) {
			idx = i;
			break;
		}
	}
	if (idx == atoms.size())
		throw runtime_error("Could not find atom " + h->name + " in residue " + rh->name + " for deletion.");
	atoms.erase(atoms.begin() + idx);
}

static void createHydronium(
	const shared_ptr<Atom>& h,
	const shared_ptr<Atom>& a,
	shared_ptr<Residue>& ra,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "HHO";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	for (size_t i = 0; i < ra->atoms.size(); ++i)
		transferCoordVel(ra->atoms[i], product->atoms[i]);

	// Last H
	product->atoms[4]->coord = h->coord;
	product->atoms[4]->velocity = h->velocity;

	hydroniumTransform(product);

	int clusterID = ra->cluster;
	product->cluster = clusterID;
	auto& cluster = clusters.at(clusterID);
	if (!cluster.water || !cluster.h3o)
		throw runtime_error("createHydronium: Cluster water, or H3O+ uninitialized.");
	*cluster.water -= 1;
	*cluster.h3o += 1;
	cluster.calculatePH();

	toRemove[ra->name].push_back(ra);
	state.residueSet[string(resType)].insert(product);
}

static void createWaterFromHydroxide(
	const shared_ptr<Atom>& h,
	const shared_ptr<Atom>& a,
	shared_ptr<Residue>& ra,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "SOL";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	// O and 1 of the H's
	transferCoordVel(ra->atoms[0], product->atoms[0]);
	transferCoordVel(ra->atoms[1], product->atoms[1]);

	// Other H, can ignore M as thats set in transform
	product->atoms[2]->coord = h->coord;
	product->atoms[2]->velocity = h->velocity;

	waterTransform(product);

	int clusterID = ra->cluster;
	product->cluster = clusterID;
	auto& cluster = clusters.at(clusterID);
	if (!cluster.water || !cluster.oh)
		throw runtime_error("createWaterFromHydroxide: Cluster water, or OH- uninitialized.");
	*cluster.water += 1;
	*cluster.oh -= 1;
	cluster.calculatePH();

	toRemove[ra->name].push_back(ra);
	state.residueSet[string(resType)].insert(product);
}

static void createAcetic(
	const shared_ptr<Atom>& h,
	const shared_ptr<Atom>& a,
	shared_ptr<Residue>& ra,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "AHX";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	for (size_t i = 0; i < ra->atoms.size(); ++i)
		transferCoordVel(ra->atoms[i], product->atoms[i]);

	if (a->name == "O1")
		swapCoordVel(product->atoms[5], product->atoms[6]);

	// Tit H
	product->atoms[7]->coord = setBondLength(product->atoms[6]->coord, h->coord, 0.1);
	product->atoms[7]->velocity = h->velocity;

	product->cluster = ra->cluster;

	toRemove[ra->name].push_back(ra);
	state.residueSet[string(resType)].insert(product);
}

static void createAmmonium(
	const shared_ptr<Atom>& h,
	const shared_ptr<Atom>& a,
	shared_ptr<Residue>& ra,
	const unordered_map<string, shared_ptr<Residue>>& templateRs,
	State& state,
	unordered_map<string, vector<shared_ptr<Residue>>>& toRemove,
	unordered_map<int, Cluster>& clusters)
{
	constexpr string_view resType = "NXH";
	auto product = cloneResidue(*templateRs.at(string(resType)));

	transferCoordVel(ra->atoms[0], product->atoms[1]);
	transferCoordVel(ra->atoms[1], product->atoms[2]);
	transferCoordVel(ra->atoms[2], product->atoms[3]);
	transferCoordVel(ra->atoms[3], product->atoms[4]);

	// Tit H
	product->atoms[0]->coord = setBondLength(product->atoms[1]->coord, h->coord, 0.1);
	product->atoms[0]->velocity = h->velocity;

	product->cluster = ra->cluster;

	toRemove[ra->name].push_back(ra);
	state.residueSet[string(resType)].insert(product);
}

static void protonateAminoAcid(
	const shared_ptr<Exchange>& e,
	const shared_ptr<Atom>& h,
	const shared_ptr<Atom>& a,
	shared_ptr<Residue>& ra)
{
	auto& atoms = ra->atoms;

	// Atom names of the newly transferred H
	// HIS has HISE/HISD so handle specially later
	static const unordered_map<string, string> newHydrogenMap = {
		{ "LYS", "HZ3" }, { "ARG", "HH12" }, { "ASP", "HD2" }, { "GLU", "HE2" },
		{ "NTERM", "H3"}, { "CTERM", "HT2"}
	};

	// Atom names in front of new H names for insertion into Residue Atom vector
	static const unordered_map<string, string> insertMap = {
		{ "LYS", "HZ2" }, { "ARG", "HH11" }, { "ASP", "OD2" }, { "GLU", "OE2" }, 
		{ "NTERM", "H2"}, { "CTERM", "OT2"}
	};

	// Actual accepting atom, for finding new H location later
	static const unordered_map<string, string> parentMap = {
		{ "LYS", "NZ" }, { "ARG", "NH1" }, { "ASP", "OD2" }, { "GLU", "OE2" },
		{ "NTERM", "N"}, { "CTERM", "OT2"}
	};

	string resType;
	string hName;
	string insertName;
	string parentName;

	if (e->isNterm) {
		resType = "NTERM";
		hName = newHydrogenMap.at(resType);
		insertName = insertMap.at(resType);
		parentName = parentMap.at(resType);
	}
	else if (e->isCterm) {
		resType = "CTERM";
		hName = newHydrogenMap.at(resType);
		insertName = insertMap.at(resType);
		parentName = parentMap.at(resType);

		if (a->name == "OT1") {
			shared_ptr<Atom> ot1, ot2;
			for (auto& atom : atoms) {
				if		(atom->name == "OT1") ot1 = a;
				else if (atom->name == "OT2") ot2 = a;
			}
			swapCoordVel(ot1, ot2);
		}
	}
	else if (ra->name == "HIS") {
		resType = "HIS";
		auto& hatoms = ra->atoms;

		if (a->name == "ND1") {
			hName = "HD1";
			insertName = a->name;
			parentName = a->name;
			swapAtoms( hatoms[7], hatoms[13]);
			swapAtoms( hatoms[8], hatoms[14]);
			swapAtoms( hatoms[9], hatoms[14]);
			swapAtoms(hatoms[10], hatoms[11]);
			swapAtoms(hatoms[11], hatoms[12]);
			swapAtoms(hatoms[12], hatoms[13]);
			swapAtoms(hatoms[13], hatoms[14]);
		}
		else if (a->name == "NE2") {
			hName = "HE2";
			insertName = a->name;
			parentName = a->name;
			swapAtoms( hatoms[7], hatoms[13]);
			swapAtoms( hatoms[8], hatoms[14]);
			swapAtoms(hatoms[10], hatoms[12]);
			swapAtoms(hatoms[11], hatoms[13]);
			swapAtoms(hatoms[12], hatoms[14]);
			swapAtoms(hatoms[13], hatoms[14]);
		}
		else
			throw runtime_error("Improper HIS accepting atom " + a->name);
	}
	else if (parentMap.contains(ra->name)) {
		resType = ra->name;
		hName = newHydrogenMap.at(resType);
		insertName = insertMap.at(resType);
		parentName = parentMap.at(resType);

		if (resType == "GLU") {
			if (a->name == "OE1") {
				shared_ptr<Atom> oe1, oe2;
				for (auto& atom : atoms) {
					if (atom->name == "OE1") oe1 = atom;
					else if (atom->name == "OE2") oe2 = atom;
				}
				swapCoordVel(oe1, oe2);
			}
		}
		else if (resType == "ASP") {
			if (a->name == "OD1") {
				shared_ptr<Atom> od1, od2;
				for (auto& atom : atoms) {
					if (atom->name == "OD1") od1 = atom;
					else if (atom->name == "OD2") od2 = atom;
				}
				swapCoordVel(od1, od2);
			}
		}
	}
	else
		throw runtime_error("protonateAminoAcid: Could not identify residue " 
			+ ra->name + " for amino acid protonation");

	// Find where to insert new hydrogen 
	auto insert_it = find_if(atoms.begin(), atoms.end(),
		[&insertName](const shared_ptr<Atom>& atom) { return atom->name == insertName; });
	if (insert_it == atoms.end())
		throw runtime_error("Could not find insertion parent " + insertName + " of residue " + ra->name + " for H insertion.");

	// Find accepting parent atom
	auto parent_it = find_if(atoms.begin(), atoms.end(),
		[&parentName](const shared_ptr<Atom>& atom) { return atom->name == parentName; });
	if (parent_it == atoms.end())
		throw runtime_error("Could not find parent " + parentName + " of residue " + ra->name + " during amino acid protonation.");

	Vec3D newCoord = setBondLength((*parent_it)->coord, h->coord, 0.1);
	atoms.insert(insert_it + 1, make_shared<Atom>(hName, 'H', newCoord, h->velocity, ra));
}

// Make template residues for inserting products
static unordered_map<string, shared_ptr<Residue>> getTemplateResidues() {
	unordered_map<string, shared_ptr<Residue>> templateRs;
	for (const auto& [resType, resGro] : Constants::templateGroNames) 
		templateRs[resType] = parseResidueGRO(resGro);
	return templateRs;
}

void doExchanges(
	const unordered_set<shared_ptr<Exchange>>& exchanges,
	State& state,
	unordered_map<int, Cluster>& clusters,
	const unordered_set<shared_ptr<Atom>>& hydrogens,
	const unordered_set<shared_ptr<Atom>>& acceptors)
{
	if (exchanges.empty()) return;

	static auto templateRs = getTemplateResidues();
	unordered_map<string, vector<shared_ptr<Residue>>> toRemove;

	for (const auto& e : exchanges) {
		auto& h = e->hydrogen;
		auto& a = e->acceptor;
		auto rh = h->parent.lock();
		auto ra = a->parent.lock();

		// Grotthuss handled specially as products == reactants
		if (e->isGrotthuss) { 
			if (rh->name == "HHO" && ra->name == "SOL") 
				doGrotthussH3O(h, a, rh, ra, state);
			else if (ra->name == "OHX" && rh->name == "SOL") 
				doGrotthussOH(h, a, rh, ra, state);
			else
				throw runtime_error("Improper Grotthuss exchange during doExchanges().");

			continue;
		}

		// Modify donor 
		if (rh->name == "SOL")
			createHydroxide(h, rh, templateRs, state, toRemove, clusters);

		else if (rh->name == "HHO")
			createWaterFromHydronium(h, rh, templateRs, state, toRemove, clusters);

		else if (rh->name == "AHX")
			createAcetate(h, rh, templateRs, state, toRemove, clusters);

		else if (rh->name == "NXH")
			createAmmonia(h, rh, templateRs, state, toRemove, clusters);

		else if (Constants::aminoAcids.contains(rh->name))
			deprotonateAminoAcid(e, h, rh);

		else
			throw runtime_error("Unknown protein donor " + rh->name + " encountered in doExchanges().");

		// Modify acceptor
		if (ra->name == "SOL")
			createHydronium(h, a, ra, templateRs, state, toRemove, clusters);

		else if (ra->name == "OHX")
			createWaterFromHydroxide(h, a, ra, templateRs, state, toRemove, clusters);

		else if (ra->name == "ATX")
			createAcetic(h, a, ra, templateRs, state, toRemove, clusters);

		else if (ra->name == "NXX")
			createAmmonium(h, a, ra, templateRs, state, toRemove, clusters);

		else if (Constants::aminoAcids.contains(ra->name))
			protonateAminoAcid(e, h, a, ra);

		else
			throw runtime_error("Unknown protein acceptor " + ra->name + " encountered in doExchanges().");
	}

	// Delete old reactants
	for (auto& [resType, residues] : toRemove) {
		auto& residueSet = state.residueSet.at(resType);
		for (auto& r : residues) {
			if (residueSet.contains(r))
				residueSet.erase(r);
			else
				throw runtime_error("Could not find Residue " + r->ID + " for deletion.");
		}
	}
}