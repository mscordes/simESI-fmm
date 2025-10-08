#include "Residue.h"
#include "Atom.h"

#include <iostream>
#include <iomanip>
#include <utility>

Residue::Residue(
    string ID_,
    string name_,
	int num_,
    char chain_)
	: ID(move(ID_)),
    name(move(name_)),
    num(num_),
	chain(chain_)
{
}

void Residue::addAtom(shared_ptr<Atom> atom) {
    atoms.emplace_back(atom);
}

void Residue::print() const
{
    cout << fixed << setprecision(2);
    cout << "\nResidue: " << ID
        << ", Name: " << name
        << ", Number: " << num
        << ", pKa: " << (pka ? to_string(*pka) : "N/A")
        << ", GPB: " << (gpb ? to_string(*gpb) : "N/A");
    if (termini) {
        cout << ", Termini pKa: " << (termini_pka ? to_string(*termini_pka) : "N/A")
            << ", Termini GPB: " << (termini_gpb ? to_string(*termini_gpb) : "N/A");
    }
    cout << ", Chain: " << chain
        << ", Cluster: " << cluster
        << "\nAtoms:\n";
    for (const auto& atom : atoms) {
        atom->print();
    }
}