#pragma once
#include "Parameters.h"
#include "Protein.h"
#include "Residue.h"
#include "Vec3D.h"

#include <memory>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>

using namespace std;

struct State {                 
    map<char, shared_ptr<Protein>> proteinMap;
    unordered_map<string, unordered_set<shared_ptr<Residue>>> residueSet;

    State() = default;

    void parseProteinPDB(const string& pdbFile);
    void parseGRO(const string& groFile);
    void update(const string& groFile);

    void writeProteinPDB(const string& fname) const;
    void writeChainGRO(const string& fname, double boxD, const Protein& protein) const;
    void writeGRO(const string& fname, const Vec3D& boxD) const;
    void writeGRO(const string& fname, double boxScalar) const;
};