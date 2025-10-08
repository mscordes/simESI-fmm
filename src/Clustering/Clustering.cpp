#include "Clustering.h"
#include "Atom.h"
#include "CoordinateManip.h"
#include "DisjointSetUnion.h"

#include <cmath>
#include <iostream>
#include <iomanip>

Cluster::Cluster(int id_, int water_)
	: id(id_), water(water_)
{
}

Cluster::Cluster(
	int id_, int water_, int h3o_, int oh_)
	: id(id_), water(water_),h3o(h3o_), oh(oh_)
{
	calculatePH();
}

void Cluster::calculatePH() {
	if (!water || !h3o || !oh)
		throw runtime_error("calculatePH: Cluster water, H3O+, or OH- uninitialized.");

	int netIons = *h3o - *oh;
	if (netIons > 0)
		ph = -log10(55.0679 * netIons / *water);
	else if (netIons < 0)
		ph = 14.0 + log10(55.0679 * abs(netIons) / *water);
	else
		ph = 7.0;
}

void Cluster::print() const {
	cout << "Cluster: " << (id ? to_string(*id) : "N/A") << ", Waters: " << (water ? to_string(*water) : "N/A") 
		<< ", H3O+: " << (h3o ? to_string(*h3o) : "N/A") << ", OH-: " << (oh ? to_string(*oh) : "N/A") 
		<< ", pH: " << setprecision(2) << (ph ? to_string(*ph) : "N/A") << '\n';
}

unordered_map<int, Cluster> clusterWaters(State& state, const WaterInfo& waterInfo) {
	vector<int> clusterIDs = computeClusters(waterInfo.Ocoords);

	const auto& waters = waterInfo.waters;
	for (size_t i = 0; i < waters.size(); ++i)
		waters[i]->cluster = clusterIDs[i];

	unordered_map<int, Cluster> clusters;
	clusters.emplace(-1, Cluster(-1, 0)); // -1 as gas-phase

	unordered_map<int, int> counts;
	for (int id : clusterIDs) 
		++counts[id];

	for (const auto& [id, count] : counts) 
		clusters.emplace(id, Cluster(id, count));
 
	return clusters;
}

// Assign cluster to hydrogens/acceptors
void clusterTitSites(
	unordered_set<shared_ptr<Atom>>& hydrogens, 
	unordered_set<shared_ptr<Atom>>& acceptors,
	State& state, 
	const WaterInfo& waterInfo, 
	double cutoff)
{
	auto it = state.residueSet.find("SOL");
	if (it == state.residueSet.end() || it->second.empty()) return;
	const auto& waters = it->second;

	for (auto& h : hydrogens) {
		auto r = h->parent.lock();

		double dist = numeric_limits<double>::infinity();
		auto o = findClosestWaterO(dist, h->coord, waterInfo);

		if (dist < cutoff) 
			r->cluster = o->parent.lock()->cluster;
		else
			r->cluster = -1; // Gas-phase
	}

	for (auto& a : acceptors) {
		auto r = a->parent.lock();

		double dist = numeric_limits<double>::infinity();
		auto h = findClosestWaterH(dist, a->coord, waterInfo);

		if (dist < cutoff)
			r->cluster = h->parent.lock()->cluster;
		else
			r->cluster = -1;
	}
}

// Create Cluster objects with info like size/pH
// Assumes that all titratable sites have already been assigned a cluster
void computeClusterPHs(unordered_map<int, Cluster>& clusters, const State& state) {	
	unordered_map<int, int> zeroes;
	for (const auto& pair : clusters) 
		zeroes[pair.first] = 0;

	unordered_map<int, int> h3o = zeroes;
	if (state.residueSet.contains("HHO")) 
		for (const auto& r : state.residueSet.at("HHO")) 
			h3o[r->cluster] += 1;

	unordered_map<int, int> oh = zeroes;
	if (state.residueSet.contains("OHX"))
		for (const auto& r : state.residueSet.at("OHX"))
			oh[r->cluster] += 1;

	// Make sure gas-phase cluster correct
	clusters[-1].h3o   = 0;
	clusters[-1].oh    = 0;
	clusters[-1].water = 0;

	for (auto& [id, cluster] : clusters) {
		cluster.h3o = h3o.at(id);
		cluster.oh = oh.at(id);
		cluster.calculatePH();
	}
}