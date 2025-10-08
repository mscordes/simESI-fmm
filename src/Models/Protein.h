#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;

struct Residue; // Forward declaration of Residue

struct Protein {
    char chain;                            
    vector<shared_ptr<Residue>> residues;

    Protein() = default;

    Protein(char chain);

	void addResidue(const shared_ptr<Residue>& residue);
	void print() const;

    double getMass() const;
};