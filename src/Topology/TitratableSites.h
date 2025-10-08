#pragma once
#include "State.h"
#include "Atom.h"

void getTitSites(
	const State& state,
	unordered_set<shared_ptr<Atom>>& hydrogens,
	unordered_set<shared_ptr<Atom>>& acceptors);