#include "FindExchanges.h"
#include "Constants.h"
#include "CoordinateManip.h"
#include "Charge.h"
#include "MathUtils.h"
#include "RNG.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <utility>

unordered_set<shared_ptr<Exchange>> findExchanges(
	const State& state,
	const unordered_set<shared_ptr<Atom>>& hydrogens,
	const unordered_set<shared_ptr<Atom>>& acceptors,
	const WaterInfo& waterInfo,
	const unordered_map<int, Cluster>& clusters,
	double temp,
	int step,
	int hop,
	double cutoff)
{
	// Get protein coordinates/charges for computing Grotthuss exchange energies
	vector<Vec3D> protCoords;
	vector<double> protCharges;
	unordered_map<string, vector<Vec3D>> sCoords;
	unordered_map<string, vector<double>> sCharges;

	const auto& residueSet = state.residueSet;
	auto h3o_it = residueSet.find("HHO");
	auto  oh_it = residueSet.find("OHX");

	if ((h3o_it != residueSet.end() && h3o_it->second.size() > 0) || 
		( oh_it != residueSet.end() &&  oh_it->second.size() > 0)) 
	{
		protCoords = getProteinCoords(state);
		protCharges = getProteinCharges(state);
		if (protCoords.size() != protCharges.size())
			throw runtime_error("Size of protien coordinates (" + to_string(protCoords.size())
				+ ") != size of charges (" + to_string(protCharges.size()) + ").");

		sCoords = getSoluteCoords(state);
		sCharges = getSoluteCharges(state);
		if (sCoords.size() != sCharges.size())
			throw runtime_error("Size of solute coordinates (" + to_string(sCoords.size())
				+ ") != size of charges (" + to_string(sCharges.size()) + ").");
	}

	unordered_set<shared_ptr<Exchange>> exchanges;
	if (!state.residueSet.contains("SOL")) return {};
	const auto& waters = state.residueSet.at("SOL");

	// For RH + H2O -> R + H3O+ exchanges
	for (auto& h : hydrogens) {
		auto rh = h->parent.lock();

		// Non-Grotthuss exchanges involving water only on hop 0
		if (hop > 0 && rh->name != "HHO") continue;

		double dist = numeric_limits<double>::infinity();
		auto o = findClosestWaterO(dist, h->coord, waterInfo); 
		if (dist < cutoff) {
			exchanges.emplace(
				make_shared<Exchange>(state, h, o, step, hop, clusters.at(rh->cluster), protCoords,
					protCharges, sCoords, sCharges, waterInfo.coords, temp));
		}
	}

	// For R + H2O -> RH + OH- exchanges
	for (auto& a : acceptors) {
		auto ra = a->parent.lock();

		// Non-Grotthuss exchanges involving water only on hop 0
		if (hop > 0 && ra->name != "OHX") continue;

		double dist = numeric_limits<double>::infinity();
		auto h = findClosestWaterH(dist, a->coord, waterInfo);
		if (dist < cutoff) {
			exchanges.emplace(
				make_shared<Exchange>(state, h, a, step, hop, clusters.at(ra->cluster), protCoords,
					protCharges, sCoords, sCharges, waterInfo.coords, temp));
		}
	}

	// All other non-water involving exchanges
	double cutoff_sq = cutoff * cutoff;
	static const auto& aminoAcids = Constants::aminoAcids;

	for (const auto& h : hydrogens) {
		const auto& hCoord = h->coord;
		for (const auto& a : acceptors) {
			if (hCoord.distance_sq(a->coord) < cutoff_sq) {

				auto rh = h->parent.lock();
				auto ra = a->parent.lock();

				if (!rh) throw runtime_error("findExchanges: Expired hydrogen residue");
				if (!ra) throw runtime_error("findExchanges: Expired acceptor residue");

				auto cit = clusters.find(rh->cluster);
				if (cit == clusters.end())
					throw runtime_error("findExchanges: Could not locate cluster " + to_string(rh->cluster) + " for residue " + rh->name);

				if (rh == ra) continue; 
				if (aminoAcids.contains(rh->name) && aminoAcids.contains(ra->name)) continue;

				exchanges.emplace(
					make_shared<Exchange>(state, h, a, step, hop, cit->second, protCoords, protCharges,
						sCoords, sCharges, waterInfo.coords, temp));
			}
		}
	}

	return exchanges;
}

void pruneExchanges(
	unordered_set<shared_ptr<Exchange>>& exchanges,
	unordered_set<shared_ptr<Residue>>& skipList, 
	RunFlags& flags,
	double temp) 
{
	if (exchanges.empty()) return;

	// MCMC sampling
	for (auto it = exchanges.begin(); it != exchanges.end(); ) {
		const auto& e = *it; 
		bool remove = (e->isGrotthuss && e->energy > 0) ||
			(!e->isGrotthuss && !MCMC(e->energy, temp));

		it = remove ? exchanges.erase(it) : next(it);
	}

	// Prevent carboxylate protonation at 180 degrees
	static unordered_set<string> carboxylates = {"ATX", "GLU", "ASP", "ARG" };
	static unordered_map<string, string> carbonNames = {
		{"ATX", "C2"}, {"GLU", "CD"}, {"ASP", "CG"}, {"ARG", "CZ"}
	};

	for (auto it = exchanges.begin(); it != exchanges.end(); ) {
		const auto& e = *it;
		auto ra = e->acceptor->parent.lock();

		bool remove = false;
		if (carboxylates.contains(ra->name)) {
			string Cname = carbonNames.at(ra->name);
			const auto& hc = e->hydrogen->coord;
			const auto& ac = e->acceptor->coord;

			for (const auto& atom : ra->atoms) {
				if (atom->name == Cname) {
					if (getAngle(atom->coord, ac, hc) < 3.0) {
						remove = true;
						break;
					}
				}
			}
		}

		it = remove ? exchanges.erase(it) : next(it);
	}

	// Prevent 'rattling' of a proton of between proton donor/acceptor pair over multiple hops
	if (!skipList.empty()) {
		for (auto it = exchanges.begin(); it != exchanges.end(); ) {
			const auto& e = *it;

			auto rh = e->hydrogen->parent.lock();
			auto ra = e->acceptor->parent.lock();

			bool remove = false;
			if (e->isGrotthuss) { // Grotthuss exchanges, only consider water
				if ((rh->name == "SOL" && skipList.contains(rh)) ||
					(ra->name == "SOL" && skipList.contains(ra)))
				{
					remove = true;
				}
			}
			else { // Non-Grotthuss exchanges
				if (skipList.contains(rh) || skipList.contains(ra))
					remove = true;
			}

			it = remove ? exchanges.erase(it) : next(it);
		}
	}

	// Find most energetically exchange per residue
	{
		unordered_map<shared_ptr<Residue>, shared_ptr<Exchange>> grotthussMin;
		unordered_map<shared_ptr<Residue>, shared_ptr<Exchange>> exchangeMin;

		auto updateMin = [](auto& m, const shared_ptr<Residue>& key, const shared_ptr<Exchange>& e) {
			auto [it, inserted] = m.try_emplace(key, e);
			if (!inserted && it->second->energy >= e->energy)
				it->second = e;
			};

		// Find minima first
		for (const auto& e : exchanges) {
			const auto& h = e->hydrogen;
			const auto& a = e->acceptor;
			auto rh = h->parent.lock();
			auto ra = a->parent.lock();

			auto& table = e->isGrotthuss ? grotthussMin : exchangeMin;
			updateMin(table, rh, e);
			updateMin(table, ra, e);
		}

		// Now remove
		for (auto it = exchanges.begin(); it != exchanges.end(); ) {
			const auto& e = *it;

			const auto& h = e->hydrogen;
			const auto& a = e->acceptor;
			auto rh = h->parent.lock();
			auto ra = a->parent.lock();

			auto& table = e->isGrotthuss ? grotthussMin : exchangeMin;

			bool keep =
				table.contains(rh) && table[rh] == e &&
				table.contains(ra) && table[ra] == e;

			it = keep ? next(it) : exchanges.erase(it);
		}

		// Now compare Grotthuss and non-Grotthuss
		for (auto it = exchanges.begin(); it != exchanges.end(); ) {
			const auto& e = *it;

			const auto& h = e->hydrogen;
			const auto& a = e->acceptor;
			auto rh = h->parent.lock();
			auto ra = a->parent.lock();

			bool keep = true;
			if (grotthussMin.contains(rh) && exchangeMin.contains(rh)) {
				if (e->isGrotthuss && exchangeMin.at(rh)->energy < 0)
					keep = false;
				else if (!e->isGrotthuss && exchangeMin.at(rh)->energy > 0)
					keep = false;
			}
			else if (grotthussMin.contains(ra) && exchangeMin.contains(ra)) {
				if (e->isGrotthuss && exchangeMin.at(ra)->energy < 0)
					keep = false;
				else if (!e->isGrotthuss && exchangeMin.at(ra)->energy > 0)
					keep = false;
			}

			it = keep ? next(it) : exchanges.erase(it);
		}
	}

	// Find most energetically favorable pH changing exchange, per water cluster
	{
		unordered_map<int, shared_ptr<Exchange>> clusterMin;
		clusterMin.reserve(exchanges.size());

		for (const auto& e : exchanges) {
			if (e->clusterID == -1) continue; // Ignore gas-phase

			auto [it, inserted] = clusterMin.try_emplace(e->clusterID, e);
			if (!inserted && e->energy < it->second->energy) {
				it->second = e;
			}
		}

		for (auto it = exchanges.begin(); it != exchanges.end(); ) {
			const auto& e = *it;
			if (e->clusterID == -1) { ++it; continue; } // Ignore gas-phase

			auto cit = clusterMin.find(e->clusterID);
			if (cit == clusterMin.end()) { ++it; continue; }

			it = (cit->second == e) ? next(it) :  exchanges.erase(it);
		}
	}

	// Update flags and skipList
	if (exchanges.size() > 0) {
		flags.exchanges = true;
		flags.createRunFile = true;

		for (const auto& e : exchanges) {
			auto rh = e->hydrogen->parent.lock();
			auto ra = e->acceptor->parent.lock();

			if (e->isGrotthuss) { // Only add water of Grotthus pair to skipList
				if (rh->name == "SOL") skipList.insert(rh);
				else if (ra->name == "SOL") skipList.insert(ra);
				else throw runtime_error("Misidentified Grotthuss pair when adding to skipList?");
			}
			else { // Non-Grotthuss
				if (Constants::aminoAcids.contains(rh->name) || Constants::aminoAcids.contains(ra->name)) {
					flags.protExchanges = true;
					if (Constants::aminoAcids.contains(rh->name)) skipList.insert(ra);
					if (Constants::aminoAcids.contains(ra->name)) skipList.insert(rh);
				}
				else {
					skipList.insert(rh);
					skipList.insert(ra);
				}
			}
		}
	}
}