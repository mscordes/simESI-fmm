#include "Charge.h"
#include "Atom.h"
#include "Constants.h"
#include "ForcefieldParser.h"
#include "StringUtils.h"
#include "TopologyUtils.h"

#include <cmath>
#include <numeric>
#include <string>
#include <fstream>
#include <unordered_set>

static unordered_map<string, double> getResCharges(
	const shared_ptr<Residue>& r,
	const unordered_map<string, unordered_map<string, double>>& chargeMap,
	bool nterm, 
	bool cterm) 
{
	static const unordered_map<string, pair<string, string>> protStateMap = {
		{"LYS", {"LYS","LSN"}}, {"ARG", {"ARG", "ARGN"}}, {"GLU", {"GLUP", "GLU"}}, {"ASP", {"ASPP", "ASP"}}
	};

	string rName;
	if (r->name == "HIS") {
			if (isAA_Protonated(r)) rName = "HSP";
			else if (isHISD(r)) rName = "HSD";
			else rName = "HSE";
	}
	else if (Constants::titAminoAcids.contains(r->name)) {
		if (isAA_Protonated(r)) rName = protStateMap.at(r->name).first;
		else					rName = protStateMap.at(r->name).second;
	}
	else rName = r->name;

	auto rIt = chargeMap.find(rName);
	if (rIt == chargeMap.end())
		throw runtime_error("Could not find FF charges for residue " + rName);
	unordered_map<string, double> rCharges = rIt->second;

	string termName;
	if (nterm) {
		if (isAA_Protonated(r, true, false)) {
			if (r->name == "GLY") termName = "GLY-NH3+";
			else if (r->name == "PRO") termName = "PRO-NH2+";
			else termName = "NH3+";
		}
		else {
			if (r->name == "GLY") termName = "GLY-NH2";
			else termName = "NH2";
		}
	}
	else if (cterm) {
		if (isAA_Protonated(r, false, true)) termName = "COOH";
		else termName = "COO-";
	}

	if (!termName.empty()) {
		auto termIt = chargeMap.find(termName);
		if (termIt == chargeMap.end())
			throw runtime_error("Could not find FF charges for terminal residue " + r->ID);
		const auto& termCharges = termIt->second;
		
		for (const auto& pair : termCharges)
			rCharges[pair.first] = pair.second;
	}

	return rCharges;
}

static void pushBackCharges(
	vector<double>& protCharges,
	const shared_ptr<Residue>& r,
	const unordered_map<string, unordered_map<string, double>>& chargeMap, 
	bool nterm, bool cterm) 
{
	const auto rCharges = getResCharges(r, chargeMap, nterm, cterm);
	for (const auto& a : r->atoms) {
		auto aIt = rCharges.find(a->name);
		if (aIt == rCharges.end())
			throw runtime_error("Could not find FF charges for atom " + a->name + " of " + r->name);
		protCharges.push_back(aIt->second);
	}
}

// Get partial charge of all protein atoms
// NOTE: State must be in proper CHARMM36 format
vector<double> getProteinCharges(const State& state) {
	static auto const chargeMap = getForceFieldCharges();
	vector<double> protCharges;

	for (const auto& [chain, p] : state.proteinMap) {

		const auto& nterm = p->residues.front();
		pushBackCharges(protCharges, nterm, chargeMap, true, false);

		for (size_t i = 1; i < p->residues.size() - 1; ++i) {
			const auto& r = p->residues[i];
			pushBackCharges(protCharges, r, chargeMap, false, false);
		}

		const auto& cterm = p->residues.back();
		pushBackCharges(protCharges, cterm, chargeMap, false, true);
	}

	return protCharges;
}

double getProteinCharge(const State& state) {
	const auto protCharges = getProteinCharges(state);
	return accumulate(protCharges.begin(), protCharges.end(), 0.0);
}

double getProteinCharge(const vector<double>& protCharges) {
	return accumulate(protCharges.begin(), protCharges.end(), 0.0);
}

// Gets partial charges of all atoms excluding protein, water, and atmosphere
unordered_map<string, vector<double>> getSoluteCharges(const State& state) {
	unordered_map<string, vector<double>> charges;
	charges.reserve(Constants::topOrder.size()); 

	for (const auto& resType : Constants::topOrder) {
		if (resType == "SOL" || resType == "NNN" || resType == "OOO") continue;

		auto rIt = state.residueSet.find(resType);
		if (rIt == state.residueSet.end()) continue;
		const auto& residues = rIt->second;
		if (residues.empty()) continue;

		auto pIt = Constants::partialsMap.find(resType);
		if (pIt == Constants::partialsMap.end())
			throw runtime_error("Could not find residue type " + resType + " in partialsMap.");
		const auto& partials = pIt->second;

		vector<double> sCharges;
		sCharges.reserve(residues.size() * partials.size()); 

		for (size_t i = 0; i < residues.size(); ++i)
			sCharges.insert(sCharges.end(), partials.begin(), partials.end());

		charges.emplace(resType, move(sCharges));
	}
	return charges;
}

double getNetCharge(const State& s) {
	auto f = [&](const string& r) {auto it = s.residueSet.find(r); return it == s.residueSet.end() ? 0 : (int)it->second.size(); };
	double protCharge = getProteinCharge(getProteinCharges(s));
	return protCharge + f("NXH") + f("HHO") - f("ATX") - f("OHX");
}

double getNetCharge(const State& s, double protCharge) {
	auto f = [&](const string& r) {auto it = s.residueSet.find(r); return it == s.residueSet.end() ? 0 : (int)it->second.size(); };
	return protCharge + f("NXH") + f("HHO") - f("ATX") - f("OHX");
}