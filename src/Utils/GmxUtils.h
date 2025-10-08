#pragma once
#include "Parameters.h"
#include "State.h"
#include <string>

void call_mdrun(const string& fname);
void modifyMDPgrps(const string& opfname, const State& state, double temperature);