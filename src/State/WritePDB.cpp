#include "State.h"
#include "Atom.h"
#include "Constants.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

// Writes a .pdb file of protein(s) in the current State
// NOTE: This only writes protein atoms, will delete any other molecules in system
void State::writeProteinPDB(const string& fname) const {
    ofstream file(fname);
    if (!file.is_open()) 
        throw runtime_error("Could not open file " + fname + " for writing.");

    for (const auto& [chain, protein] : proteinMap) {
        int resCount = 0;
        int atomCount = 0;

        for (const auto& residuePtr : protein->residues) {
            Residue& residue = *residuePtr;

            if (Constants::aminoAcids.find(residue.name) == Constants::aminoAcids.end()) continue;
            ++resCount;

            for (const auto& atomPtr : residue.atoms) {
                ++atomCount;
                Atom& atom = *atomPtr;

                file << "ATOM" << right << setw(7) << atomCount
                    << "  " << left << setw(4) << atom.name
                    << left << setw(4) << residue.name
                    << chain << right << setw(4) << resCount
                    << fixed << setprecision(3) << right
                    << setw(12) << atom.coord.x * 10
                    << setw(8) << atom.coord.y * 10
                    << setw(8) << atom.coord.z * 10
                    << setw(6) << "1.00" << setw(6) << "0.00"
                    << setw(12) << atom.element << '\n';
            }
        }
        file << "TER\n";
    }
    file << "ENDMDL\n";
    file.close();
}