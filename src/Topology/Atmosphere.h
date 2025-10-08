#pragma once
#include "State.h"

int getNumGas(const double boxD);
unordered_map<string, int> getGasNums();
void seedAtmosphere(State& state);