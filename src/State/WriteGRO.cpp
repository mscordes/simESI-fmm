#include "State.h"
#include "Atom.h"
#include "Constants.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

static inline void writeGROatom(
    ofstream& file,
    const shared_ptr<Atom>& atom, 
    const string& resName,
    int& atomCount, 
    int& resCount)
{
    if (++atomCount > 99999) atomCount = 1;
    if (resCount > 99999) resCount = 1;

    char buffer[256]; // plenty of space for one line
    int len = snprintf(buffer, sizeof(buffer),
        "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
        resCount,
        resName.c_str(),
        atom->name.c_str(),
        atomCount,
        atom->coord.x, atom->coord.y, atom->coord.z,
        atom->velocity.x, atom->velocity.y, atom->velocity.z);

    file.write(buffer, len);
}

// Write a .gro of a given chain from the current State
void State::writeChainGRO(const string& fname, double boxD, const Protein& protein) const {
    ofstream file(fname);
    if (!file.is_open()) 
        throw runtime_error("Could not open file " + fname + " for writing.");

    file << "simESI-fmm generated .gro file\n";
    streampos countPos = file.tellp();
    file << setw(10) << "\n";

    int resCount = 0;
    int trueAtomCount = 0;
    int atomCount = 0;

    for (const auto& residuePtr : protein.residues) {
        ++resCount;
        const shared_ptr<Residue>& residue = residuePtr;
        for (const auto& atom : residue->atoms) {
            ++trueAtomCount;
            writeGROatom(file, atom, residue->name, atomCount, resCount);
        }
    }

    file << fixed << setprecision(5) << "  " << boxD << "  " << boxD << "  " << boxD << '\n';
    file.seekp(countPos);
    file << trueAtomCount;
    file.close();
}

// Write a complete .gro from a given State
void State::writeGRO(const string& fname, const Vec3D& boxD) const {
    ofstream file(fname);
    if (!file.is_open()) 
        throw runtime_error("Could not open file " + fname + " for writing.");

    file << "simESI-fmm generated .gro file\n";
    streampos countPos = file.tellp();
    file << setw(10) << "\n";

    int resCount = 0;
    int trueAtomCount = 0;
    int atomCount = 0;

    for (const auto& [chain, protein] : proteinMap) { // Write protein atoms first
        for (const auto& residue: protein->residues) {
            ++resCount;
            for (const auto& atom : residue->atoms) {
                ++trueAtomCount;
                writeGROatom(file, atom, residue->name, atomCount, resCount);
            }
        }
    }

    for (const string& resType : Constants::topOrder) { // Non-protein atoms
        auto it = residueSet.find(resType);
        if (it == residueSet.end()) continue;
        const auto& residues = it->second;
        for (const auto& residue : residues) {
            ++resCount;
            for (const auto& atom : residue->atoms) {
                ++trueAtomCount;
                writeGROatom(file, atom, residue->name, atomCount, resCount);
            }
        }
    }

    file << fixed << setprecision(5) << "  " << boxD.x << "  " << boxD.y << "  " << boxD.z << '\n';
    file.seekp(countPos);
    file << trueAtomCount;
    file.close();
}

void State::writeGRO(const string& fname, double boxScalar) const {
    writeGRO(fname, Vec3D{ boxScalar, boxScalar, boxScalar });
}