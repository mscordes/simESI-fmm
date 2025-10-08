#pragma once
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

using namespace std;

struct Atom; // Forward declaration of Atom

struct Residue {
    string ID;
    string name;
    int    num;
    bool   termini{false}; // Is N or C-termini
    optional<double> pka;
    optional<double> termini_pka; // Must differentiate termini pKa if termini titratable (ie an Nter LYS)
    optional<double> gpb;
    optional<double> termini_gpb;
    int    cluster{-1}; // -1 as default in gas-phase 
    char   chain;
    vector<shared_ptr<Atom>> atoms;

    Residue() = default;

    Residue(string ID,
            string name,
            int num,
            char chain);

    void addAtom(shared_ptr<Atom> atom);
	void print() const;
};