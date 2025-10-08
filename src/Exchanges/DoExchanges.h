#pragma once
#include "State.h"
#include "Exchange.h"

void doExchanges(
	const unordered_set<shared_ptr<Exchange>>& exchanges,
	State& state,
	unordered_map<int, Cluster>& clusters,
	const unordered_set<shared_ptr<Atom>>& hydrogens,
	const unordered_set<shared_ptr<Atom>>& acceptors);