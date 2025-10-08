#pragma once
#include "State.h"
#include "System.h"

unordered_set<shared_ptr<Residue>> removeEvaporated(State& state, RunFlags& flags, double cutoff = 5.0);