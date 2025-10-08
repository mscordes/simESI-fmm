#include "Evaporation.h"
#include "Atom.h"
#include "Atmosphere.h"
#include "Constants.h"
#include "CoordinateManip.h"
#include "RNG.h"

#include <cstdlib>
#include <limits>

static unordered_map<string, shared_ptr<Residue>> getTemplateGasResidues() {
    const auto gases = getGasNums();
    unordered_map<string, shared_ptr<Residue>> templateRs;
    for (const auto& [gasType, _] : gases) {
        const auto& gasGro = Constants::templateGroNames.at(gasType);
        templateRs[gasType] = parseResidueGRO(gasGro);
    }
    return templateRs;
}

static unordered_set<shared_ptr<Residue>> maintainPressure(
    State& state,
    RunFlags& flags)
{
    if (!Parameters::Get().isHumid()) return {};

    static const auto gases = getGasNums();
    static const auto templateGases = getTemplateGasResidues();

    const auto skeleton = getProteinSkeletonCoords(state);
    auto allCoords = getAllCoords(state);
    auto gWaters = getGaseousWaters(state);

    for (const auto& [gasType, n] : gases) {
        auto it = state.residueSet.find(gasType);
        if (it == state.residueSet.end()) continue;

        auto& gas = it->second;
        if (gas.empty()) continue;
        if (gasType == "SOL" && gWaters.empty()) continue;

        int delta;
        if (gasType == "SOL") delta = n - gWaters.size();
        else                  delta = n - gas.size();

        if (delta == 0) continue;

        // +/- 5% slop
        int slop = (n == 0) ? 0 : static_cast<int>(round(n * 0.05));
        if (abs(delta) < slop) continue;

        flags.createRunFile = true;
        flags.exchanges = true;

        // Need to delete gaseous molecs
        if (delta < 0) {
            delta = abs(delta);

            vector<shared_ptr<Residue>> gasVec;
            if (gasType == "SOL") {
                gasVec.reserve(gWaters.size());
                for (const auto& r : gWaters)
                    gasVec.push_back(r);
            }
            else {
                gasVec.reserve(gas.size());
                for (const auto& r : gas)
                    gasVec.push_back(r);
            }

            unordered_set<shared_ptr<Residue>> toDelete;
            toDelete.reserve(delta);
            uniform_int_distribution<size_t> d(0, gasVec.size() - 1);

            for (int i = 0; i < delta; ++i)
                toDelete.insert(gasVec[d(RNG::gen)]);
            
            for (const auto& d : toDelete) {
                gas.erase(d);
                if (gasType == "SOL") gWaters.erase(d);
            }
        }

        // Need to seed new gas molecules
        else if (delta > 0) {
            static uniform_real_distribution<double> ld(0.5, Parameters::Get().getBoxSize() - 0.5f);

            for (int i = 0; i < delta; ++i) {
                auto r = cloneResidue(*templateGases.at(string(gasType)));
                int attempts = 0;

                while (true) {
                    if (++attempts > delta * 100)
                        throw runtime_error("maintainPressure: Could not find suitable locations to seed new gas molecules");

                    Vec3D loc = { ld(RNG::gen), ld(RNG::gen), ld(RNG::gen) };
                    for (auto& a : r->atoms)
                        a->coord = a->coord + loc;

                    bool inGas = true;
                    const auto& fc = r->atoms[0]->coord;
                    for (const auto& c : allCoords) {
                        if (fc.distance_sq(c) < 1.0) {
                            for (const auto& a : r->atoms) {
                                if (a->coord.distance_sq(c) < 0.25) { // 0.5 nm cutoff
                                    inGas = false;
                                    goto endloop;
                                }
                            }
                        }
                    }
                    endloop:

                    if (inGas) {
                        gas.insert(r);
                        if (gasType == "SOL") 
                            gWaters.insert(r);
                        for (const auto& a : r->atoms)
                            allCoords.push_back(a->coord);
                        break;
                    }
                }
            }
        }
    }

    return gWaters;
}

unordered_set<shared_ptr<Residue>> removeEvaporated(State& state, RunFlags& flags, double cutoff) {
    static const auto gases = getGasNums();
    vector<Vec3D> skeleton = getProteinSkeletonCoords(state);
    const double cutoff_sq = cutoff * cutoff;
    bool deletions = false;

    unordered_set<shared_ptr<Residue>>gWaters = maintainPressure(state, flags);

    for (const auto& resType : Constants::topOrder) {
        if (gases.contains(resType)) continue;

        auto it = state.residueSet.find(resType);
        if (it == state.residueSet.end()) continue; 
  
        unordered_set<shared_ptr<Residue>> toDelete;
        for (const auto& r : it->second) {
            const auto& rc = r->atoms[0]->coord;

            double min_d2 = numeric_limits<double>::max();
            for (const auto& bc : skeleton) {
                double d2 = rc.distance_sq(bc);
                if (d2 < min_d2) {
                    min_d2 = d2;
                    if (d2 < cutoff_sq)
                        break;
                }
            }

            if (min_d2 > cutoff_sq) {
                toDelete.insert(r);
                deletions = true;
            }
        }

        for (const auto& d : toDelete)
            state.residueSet.at(resType).erase(d);
    }

    if (deletions) {
        flags.createRunFile = true;
        flags.exchanges = true;
    }

    return gWaters;
}