#pragma once
#include "State.h"

#include <map>
#include <vector>
#include <string>

using namespace std;

pair<string, string> getNtermPdb2GmxCodes(const shared_ptr<Residue>& residue);
bool isAA_Protonated(const shared_ptr<Residue>& residue, bool isNterm = false, bool isCterm = false);
bool isHISD(const shared_ptr<Residue>& his);
map<char, vector<string>> setPdb2GmxInputs(State& state, double ph);
void createTOP(const string& fname, const State& state);
map<char, vector<string>> getPdb2GmxInputs(State& state);
void call_pdb2gmx(const string& topName, const map<char, vector<string>>& gmxInputs, State& state, bool keepGro = false);
void testTOP(const State& s, string top);