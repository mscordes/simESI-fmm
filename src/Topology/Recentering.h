#pragma once
#include "State.h"
#include "System.h"

void recenter(
	State& state,
	unordered_set<shared_ptr<Residue>>& gWaters,
	RunFlags& flags, 
	double cutoff = 5.0);