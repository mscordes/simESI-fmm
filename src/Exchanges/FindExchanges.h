#pragma once
#include "Atom.h"
#include "Clustering.h"
#include "Exchange.h"
#include "State.h"

unordered_set<shared_ptr<Exchange>> findExchanges(
	const State& state,
	const unordered_set<shared_ptr<Atom>>& hydrogens,
	const unordered_set<shared_ptr<Atom>>& acceptors,
	const WaterInfo& waterInfo,
	const unordered_map<int, Cluster>& clusters,
	double temp, int step, int hop,
	double cutoff = 0.25);

void pruneExchanges(
	unordered_set<shared_ptr<Exchange>>& exchanges,
	unordered_set<shared_ptr<Residue>>& skipList, 
	RunFlags& flags,
	double temp);