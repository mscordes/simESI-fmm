#include "Exchange.h"
#include "Constants.h"
#include "CoordinateManip.h"
#include "RNG.h"

#include <cmath>
#include <iostream>
#include <iomanip>

Exchange::Exchange(
	const State& state,
	const shared_ptr<Atom>& hydrogen_,
	const shared_ptr<Atom>& acceptor_,
	int step_, int hop_,
	const Cluster& cluster,
	const vector<Vec3D>& protCoords,
	const vector<double>& protCharges,
	const unordered_map<string, vector<Vec3D>>& sCoords,
	const unordered_map<string, vector<double>>& sCharges,
	const vector<WaterQuad>& waterCoords,
	double temp_)
	: hydrogen(hydrogen_), acceptor(acceptor_),
	step(step_), hop(hop_), temp(temp_)
{
	if (!cluster.id || !cluster.water || !cluster.ph)
		throw runtime_error("Exchange::Exchange: Cluster ID, water, or pH uninitialized.");
	clusterID = *cluster.id;
	nearWaters = *cluster.water;
	ph = *cluster.ph;

	auto rh = hydrogen_->parent.lock();
	auto ra = acceptor_->parent.lock();

	try {
		// --- Donor ---
		if (rh->name == "SOL") {
			--nearWaters; // Don't count reactant water in count
			h_pka = ph; // water pKa ~ cluster pH
			h_gpb = Constants::GPBs.at("OHX");
		}
		else if (rh->termini &&
			(hydrogen_->name == "HT2" || Constants::hydrogensMap.at("NTERM").contains(hydrogen_->name)))
		{
			if (Constants::hydrogensMap.at("NTERM").contains(hydrogen_->name)) isNterm = true;
			else if (hydrogen_->name == "HT2") isCterm = true;
			else throw runtime_error("Misidentified donor termini during Exchange construction.");

			h_pka = rh->termini_pka;
			h_gpb = rh->termini_gpb;
		}
		else {
			h_pka = rh->pka;
			h_gpb = rh->gpb;
		}

		// --- Acceptor ---
		if (ra->name == "SOL") {
			--nearWaters;
			a_pka = ph;
			a_gpb = Constants::GPBs.at("HHO");
		}
		else if (ra->termini && 
			(acceptor_->name == "N" || Constants::acceptorsMap.at("CTERM").contains(acceptor_->name))) 
		{
			if (acceptor_->name == "N") isNterm = true;
			else if (Constants::acceptorsMap.at("CTERM").contains(acceptor_->name)) isCterm = true;
			else throw runtime_error("Misidentified acceptor termini during Exchange construction.");

			a_pka = ra->termini_pka;
			a_gpb = ra->termini_gpb;
		}
		else {
			a_pka = ra->pka;
			a_gpb = ra->gpb;
		}
	}
	catch (const out_of_range& e) {
		throw runtime_error(string("GPB lookup failed: ") + e.what());
	}

	if (!h_pka || !h_gpb)
		throw runtime_error("calculateEnergy: pKa or GPB not initialized for " + rh->name + ".");
	if (!a_pka || !a_gpb)
		throw runtime_error("calculateEnergy: pKa or GPB not initialized for " + ra->name + ".");

	isGrotthuss = (rh->name == "HHO" && ra->name == "SOL") || (ra->name == "OHX" && rh->name == "SOL");
	if (isGrotthuss) 
		energy = calcGrotthussEnergy(state, protCoords, protCharges, sCoords, sCharges, waterCoords);
	else 
		energy = calcExchangeEnergy();
}

// For intramolecular (mobile proton) protein transfers
Exchange::Exchange(
	const State& state,
	const shared_ptr<Atom>& hydrogen_,
	const shared_ptr<Atom>& acceptor_,
	double energy_,
	int step_, 
	double temp_)
	: hydrogen(hydrogen_), acceptor(acceptor_), energy(energy_), 
	step(step_), temp(temp_)
{
	isGrotthuss = false; // These don't matter for these types of exchanges
	hop = 0; 

	clusterID = -1; // Gas-phase only
	nearWaters = 0;
	ph = 7.0;

	auto rh = hydrogen_->parent.lock();
	auto ra = acceptor_->parent.lock();

	try {
		// --- Donor ---
		if (rh->termini &&
			(hydrogen_->name == "HT2" || Constants::hydrogensMap.at("NTERM").contains(hydrogen_->name)))
		{
			if (Constants::hydrogensMap.at("NTERM").contains(hydrogen_->name)) isNterm = true;
			else if (hydrogen_->name == "HT2") isCterm = true;
			else throw runtime_error("Misidentified donor termini during Exchange construction.");

			h_pka = rh->termini_pka;
			h_gpb = rh->termini_gpb;
		}
		else {
			h_pka = rh->pka;
			h_gpb = rh->gpb;
		}

		// --- Acceptor ---
		if (ra->termini &&
			(acceptor_->name == "N" || Constants::acceptorsMap.at("CTERM").contains(acceptor_->name)))
		{
			if (acceptor_->name == "N") isNterm = true;
			else if (Constants::acceptorsMap.at("CTERM").contains(acceptor_->name)) isCterm = true;
			else throw runtime_error("Misidentified acceptor termini during Exchange construction.");

			a_pka = ra->termini_pka;
			a_gpb = ra->termini_gpb;
		}
		else {
			a_pka = ra->pka;
			a_gpb = ra->gpb;
		}
	}
	catch (const out_of_range& e) {
		throw runtime_error(string("GPB lookup failed: ") + e.what());
	}

	if (!h_pka || !h_gpb)
		throw runtime_error("calculateEnergy: pKa or GPB not initialized for " + rh->name + ".");
	if (!a_pka || !a_gpb)
		throw runtime_error("calculateEnergy: pKa or GPB not initialized for " + ra->name + ".");
}

void Exchange::print(ostream& os) const {
	if (!hydrogen || !acceptor) 
		throw runtime_error("[Exchange] Invalid: missing hydrogen or acceptor pointer.\n");

	const auto& h = hydrogen;
	const auto& a = acceptor;

	auto rh = h->parent.lock();
	auto ra = a->parent.lock();

	if (!rh || !ra) {
		cerr << "[Exchange] Invalid: parent residue expired for "
			<< (rh ? "" : "hydrogen") << (ra ? "" : " and acceptor") << ".\n";
		throw runtime_error("[Exchange] Invalid\n");
	}

	os << "Step: " << setw(6) << step
		<< ", Hop: " << setw(2) << hop
		<< ", E: " << setw(8) << fixed << setprecision(3) << energy
		<< ", Cluster: " << setw(6) << clusterID
		<< ", Waters: " << setw(6) << nearWaters
		<< ", pH: " << setw(4) << setprecision(2) << ph;

	os << "\n    " << setw(5) << h->name
		<< setw(10) << rh->ID
		<< ", pKa: " << setw(5) << fixed << setprecision(2) << *h_pka
		<< ", GPB: " << setw(7) << setprecision(2) << *h_gpb
		<< ", [" << setw(7) << fixed << setprecision(3) << h->coord.x
		<< ", " << setw(7) << h->coord.y
		<< ", " << setw(7) << h->coord.z << "]\n";

	os << "    " << setw(5) << a->name
		<< setw(10) << ra->ID
		<< ", pKa: " << setw(5) << fixed << setprecision(2) << *a_pka
		<< ", GPB: " << setw(7) << setprecision(2) << *a_gpb
		<< ", [" << setw(7) << fixed << setprecision(3) << a->coord.x
		<< ", " << setw(7) << a->coord.y
		<< ", " << setw(7) << a->coord.z << "]\n\n";
}

double Exchange::calcExchangeEnergy() const {
	constexpr double kSol = 0.019144; // ln(10)*kb
	constexpr double kExp = 0.30312; // Fit from Kumar, Phys. Chem. Chem. Phys. (2022) 24,30, 18236-18244

	auto gp_correct = [&](double gpb, double pka) { // Gas-phase correction
		if (nearWaters < 1) return gpb;
		double sol = kSol * pka * temp;
		return ((gpb - sol) * exp(-kExp * nearWaters)) + sol;
		};

	double energy;
	if (nearWaters > 30)
		energy = kSol * (*h_pka - *a_pka) * temp;
	else if (nearWaters < 1)
		energy = *h_gpb - *a_gpb;
	else
		energy = gp_correct(*h_gpb, *h_pka) - gp_correct(*a_gpb, *a_pka);

	static uniform_real_distribution<double> d(-0.01, 0.01); // For comparisons
	energy += d(RNG::gen);

	if (isinf(energy)) {
		cerr << "inf energy " + to_string(energy) + " encountered during exchange energy calculation " +
			"for exchange, \n";
		Exchange::print(cout);
		throw runtime_error("calcExchangeEnergy: inf energy exchange");
	}

	return energy;
}

// Couloumb potential at a given location
static double couloumbEnergy(
	const Vec3D& initLoc, 
	const Vec3D& finalLoc, 
	int charge,
	const vector<Vec3D>& protCoords,
	const vector<double>& protCharges,
	const unordered_map<string, vector<Vec3D>>& sCoords,
	const unordered_map<string, vector<double>>& sCharges,
	const vector<WaterQuad>& waterCoords, 
	const shared_ptr<Residue>& rh, 
	const shared_ptr<Residue>& ra, 
	const State& state)
{
	double initEnergy = 0.0;
	double finalEnergy = 0.0;

	// Water contribution
	static double Mcharge = Constants::partialsMap.at("SOL")[3]; // Virtual site
	static double Hcharge = Constants::partialsMap.at("SOL")[1];
	static double cutoff = 0.50; // 0.50 nm cutoff

	for (const auto& w : waterCoords) {

		// Avoid waters within cutoff of reference point to broadly sample ES environment
		if (initLoc.distance(w.O) < cutoff) continue;

		initEnergy += Mcharge / initLoc.distance(w.M);
		initEnergy += Hcharge / initLoc.distance(w.H1);
		initEnergy += Hcharge / initLoc.distance(w.H2);

		finalEnergy += Mcharge / finalLoc.distance(w.M);
		finalEnergy += Hcharge / finalLoc.distance(w.H1);
		finalEnergy += Hcharge / finalLoc.distance(w.H2);
	}

	// Protein Contribution
	for (size_t i = 0; i < protCoords.size(); ++i) {
		initEnergy  += protCharges[i] /  initLoc.distance(protCoords[i]);
		finalEnergy += protCharges[i] / finalLoc.distance(protCoords[i]);
	}

	// Solute contribution, must avoid the H3O+/OH- involved in the exchange
	for (const auto& [resType, coords] : sCoords) {
		auto pIt = sCharges.find(resType);
		if (pIt == sCharges.end())
			throw runtime_error("couloumbEnergy: Could not find residue type " + resType + " in sCharges.");
		const auto& charges = pIt->second;

		// Get solute coords ignoring Grotthuss H3O+/OH- if necessary
		if ((resType == "OHX" && ra->name == "OHX") || (resType == "HHO" && rh->name == "HHO")) {
			auto it = state.residueSet.find(resType);
			if (it == state.residueSet.end())
				throw runtime_error("couloumbEnergy: Could not find residue type " + resType + " in residueSet.");
			const auto& residues = it->second;

			int idx;
			if (resType == "OHX" && ra->name == "OHX") idx = 0;
			else if (resType == "HHO" && rh->name == "HHO") idx = 3;
			else throw runtime_error("couloumbEnergy: Bad exchange type during solute ES.");

			for (const auto& r : residues) {
				const auto& atoms = r->atoms;
				double rDist = initLoc.distance(atoms[idx]->coord);
				if (rDist < 0.05) continue;
				
				for (size_t i = 0; i < atoms.size(); ++i) {
					// Charges just repeating pattern so indexing here is fine 
					initEnergy  += charges[i] /  initLoc.distance(atoms[i]->coord);
					finalEnergy += charges[i] / finalLoc.distance(atoms[i]->coord);
				}
			}
		}
		else {
			for (size_t i = 0; i < coords.size(); ++i) {
				initEnergy  += charges[i] /  initLoc.distance(coords[i]);
				finalEnergy += charges[i] / finalLoc.distance(coords[i]);
			}
		}
	}

	// Finish ES calc and return delta
	static const double ke = Constants::ke;
	initEnergy  *= charge * ke;
	finalEnergy *= charge * ke; 
	return finalEnergy - initEnergy; // kJ/mol
}

// Compute delta Couloumb energy for Grotthuss mechanism 
double Exchange::calcGrotthussEnergy(
	const State& state, 
	const vector<Vec3D>& protCoords,
	const vector<double>& protCharges,
	const unordered_map<string, vector<Vec3D>>& sCoords,
	const unordered_map<string, vector<double>>& sCharges,
	const vector<WaterQuad>& waterCoords) const 
{
	const auto& h = hydrogen;
	const auto& a = acceptor;
	auto rh = h->parent.lock();
	auto ra = a->parent.lock();

	double energy;
	if (rh->name == "HHO" && ra->name == "SOL") { // H3O+ + H2O -> H2O + H3O+
		const Vec3D& initLoc = rh->atoms[3]->coord;	// Virtual site of H3O+
		const Vec3D& finalLoc = ra->atoms[3]->coord;	// Virtual site of H2O
		int charge = 1; // +1 charge of H3O+

		energy = couloumbEnergy(initLoc, finalLoc, charge, protCoords, protCharges, 
			sCoords, sCharges, waterCoords, rh, ra, state);
		energy -= 25.0; // Diffusion coefficient correction
	}
	else if (ra->name == "OHX" && rh->name == "SOL") { // OH- + H2O -> H2O + OH-
		const Vec3D& initLoc = a->coord;				// O of OH-
		const Vec3D& finalLoc = rh->atoms[3]->coord;	// Virtual site of H2O
		int charge = -1; // -1 charge of OH-

		energy = couloumbEnergy(initLoc, finalLoc, charge, protCoords, protCharges, 
			sCoords, sCharges, waterCoords, rh, ra, state);
		if (hop > 2) energy += 9999.0; // OH- hops less frequently than H3O+
	}
	else throw runtime_error("Invalid Grotthuss pair type when computing exchange energy.");

	if (isinf(energy)) {
		cerr << "inf energy " + to_string(energy) + " encountered during Grotthuss exchange energy calculation.\n";
		Exchange::print(cout);
		throw runtime_error("calcGrotthussEnergy: inf energy exchange");
	}

	return energy;
}