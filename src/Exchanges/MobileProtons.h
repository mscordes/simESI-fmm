#pragma once
#include "Atom.h"
#include "Exchange.h"
#include "ExchangeWriter.h"
#include "State.h"

struct MobileProtonSite {
    shared_ptr<Residue> residue;
    Vec3D               coord;
    bool                protonated;
    double              charge;
    int                 idx;
    bool                nterm;
    bool                cterm;

    MobileProtonSite(
        const shared_ptr<Residue>& residue,
        int idx, bool nterm, bool cterm);
};

void mobileProtons(
    State& state,
    RunFlags& flags,
    ExchangeWriter& exchangeWriter,
    const unordered_map<int, Cluster>& clusters,
    const int step,
    double temp);