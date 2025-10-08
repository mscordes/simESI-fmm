#include "Protein.h"
#include "Atom.h"
#include "Constants.h"
#include "Residue.h"

Protein::Protein(char chain_) : chain(chain_) {}

void Protein::addResidue(const shared_ptr<Residue>& residue) {
    residues.push_back(residue);
}

void Protein::print() const {
    cout << "Protein: " << chain;
    if (!residues.empty()) {
        cout << ", from residue " << residues.front()->ID << " to " << residues.back()->ID << ".\n";
    }
}

double Protein::getMass() const {
    double mass = 0;
    for (const auto& residue : residues) {
        for (const auto& atom : residue->atoms) {
            const auto it = Constants::atomMasses.find(atom->element);
            if (it == Constants::atomMasses.end()) 
                throw runtime_error("Unknown atom element: " + atom->element);
            mass += it->second;
        }
    }
    return mass;
}