#pragma once
#include "State.h"
#include "CoordinateManip.h"

#include <optional>

struct Cluster {
	optional<int>	 id;
	optional<int>    water;
	optional<int>    h3o;
	optional<int>    oh;
	optional<double> ph;

	Cluster() = default;
	Cluster(int id, int water);
	Cluster(int id, int water, int h3o, int oh);
	void calculatePH();
	void print() const;
};

unordered_map<int, Cluster> clusterWaters(State& state, const WaterInfo& waterInfo);
void clusterTitSites(
	unordered_set<shared_ptr<Atom>>& hydrogens,
	unordered_set<shared_ptr<Atom>>& acceptors,
	State& state,
	const WaterInfo& waterInfo, 
	double cutoff = 0.40);
void computeClusterPHs(unordered_map<int, Cluster>& clusters, const State& state);