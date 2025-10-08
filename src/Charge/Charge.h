#pragma once
#include "State.h"
#include <vector>

using namespace std;

vector<double> getProteinCharges(const State& state);
double getProteinCharge(const State& state);
double getProteinCharge(const vector<double>& protCharges);
unordered_map<string, vector<double>> getSoluteCharges(const State& state);
double getNetCharge(const State& s);
double getNetCharge(const State& s, double protCharge);