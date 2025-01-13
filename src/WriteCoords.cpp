#include "Core.h"
#include <fstream>
#include <iomanip>
#include <memory>

namespace Core {

    // Writes a .pdb file of a given protein(s) 
    // NOTE: This only writes protein atoms, will delete any other molecules in system
    void writePDB(const std::string& filename, std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms,
        const std::array<float, 3>& box_vectors) {
        
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::exit(1);
        }

        // Write box vectors
        file << "CRYST1  "
            << std::fixed << std::setprecision(2)
            << std::setw(9) << box_vectors[0] * 10
            << std::setw(9) << box_vectors[1] * 10
            << std::setw(9) << box_vectors[2] * 10
            << "  90.00  90.00  90.00 P 1           1\n";

        file << "MODEL        1\n";

        const std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        // Iterate over proteins
        int protCount = 0;
        for (const auto& pair : proteinAtoms) {
            const std::string& chain = pair.first;

            // Iterate over atoms in the protein
            int atomCount = 1;
            int resCount = 0;
            std::string temp_res_id = "0";
            for (const auto& atom : pair.second) {

                // Determine if new residue for residue numbering
                if (atom->res_id != temp_res_id) {
                    resCount += 1;
                    temp_res_id = atom->res_id;
                }

                // Start building the ATOM line
                file << "ATOM" << std::right << std::setw(7) << atomCount;  // Atom number

                // Atom name
                if (atom->atom_name.size() > 3) {
                    file << " " << std::left << std::setw(4) << atom->atom_name << " ";
                }
                else {
                    file << "  " << std::left << std::setw(4) << atom->atom_name;  // Atom name
                }

                file << std::left << std::setw(4) << atom->res_name;  // Residue name

                file << alphabet[protCount];  // Chain 

                file << std::right << std::setw(4) << resCount;  // Residue number

                // Write coordinates
                file << std::fixed << std::setprecision(3)
                    << std::right << std::setw(12) << atom->coord[0] * 10
                    << std::setw(8) << atom->coord[1] * 10
                    << std::setw(8) << atom->coord[2] * 10;

                // Write occupancy and temp factors
                file << std::setw(6) << "1.00"; 
                file << std::setw(6) << "0.00"; 

                file << std::setw(12) << atom->element;  // Element symbol

                file << '\n';  // Newline after each ATOM entry
                atomCount++;
            }

            file << "TER\n";  // Terminator after each protein's atoms
            protCount += 1;
        }

        file << "ENDMDL\n";  // End of model

        file.close(); 
    }

    // Writes a .gro file given an input list of Atom objects 
    void writeGRO(const std::string& filename, const std::vector<Atom>& atoms, const std::array<float, 3>& box_vectors) {
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::exit(1);
        }

        // Write title and number of atoms
        file << "simESI generated .gro file\n" << std::to_string(atoms.size()) << "\n";

        // Iterate over atoms
        int atomCount = 1;
        int resCount = 0;
        std::string temp_res_id = "0";
        for (const auto& atom : atoms) {
            // Determine if new residue for residue numbering
            if (atom.res_id != temp_res_id) {
                resCount += 1;

                // Correct for formatting overflow at 100,000 residues
                if (resCount == 100000) {
                    resCount = 1;
                }

                temp_res_id = atom.res_id;
            }

            // Start building each atom line
            file << std::setw(5) << resCount  // Residue number 
                << std::left << std::setw(5) << atom.res_name  // Residue name
                << std::right << std::setw(5) << atom.atom_name  // Atom name
                << std::setw(5) << atomCount;  // Atom number

            // Write coordinates
            file << std::fixed << std::setprecision(3)
                << std::setw(8) << atom.coord[0]
                << std::setw(8) << atom.coord[1]
                << std::setw(8) << atom.coord[2];

            // Write velocities
            file << std::fixed << std::setprecision(4)
                << std::setw(8) << atom.velocity[0]
                << std::setw(8) << atom.velocity[1]
                << std::setw(8) << atom.velocity[2];

            file << '\n';  // Newline after each ATOM entry

            // Count up iterator
            ++atomCount;

            // Correct for formatting overflow at 100,000 atoms
			if (atomCount == 100000) {
				atomCount = 1;
			}
        }
        
        // Box vectors
		file << std::fixed << std::setprecision(5)
			<< "  " << box_vectors[0] << "  " << box_vectors[1] << "  " << box_vectors[2] << '\n';

        file.close(); 
    }

    // Overload, writes a .gro file given an input list of pointers to Atom objects 
    void writeGRO(const std::string& filename, const std::vector<std::shared_ptr<Atom>>& atoms, const std::array<float, 3>& box_vectors) {
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::exit(1);
        }

        // Write title and number of atoms
        file << "simESI generated .gro file\n" << std::to_string(atoms.size()) << "\n";

        // Iterate over atoms
        int atomCount = 1;
        int resCount = 0;
        std::string temp_res_id = "0";
        for (const auto& atom : atoms) {
            // Determine if new residue for residue numbering
            if (atom->res_id != temp_res_id) {
                resCount += 1;

                // Correct for formatting overflow at 100,000 residues
                if (resCount == 100000) {
                    resCount = 1;
                }

                temp_res_id = atom->res_id;
            }

            // Start building each atom line
            file << std::setw(5) << resCount  // Residue number 
                << std::left << std::setw(5) << atom->res_name  // Residue name
                << std::right << std::setw(5) << atom->atom_name  // Atom name
                << std::setw(5) << atomCount;  // Atom number

            // Write coordinates
            file << std::fixed << std::setprecision(3)
                << std::setw(8) << atom->coord[0]
                << std::setw(8) << atom->coord[1]
                << std::setw(8) << atom->coord[2];

            // Write velocities
            file << std::fixed << std::setprecision(4)
                << std::setw(8) << atom->velocity[0]
                << std::setw(8) << atom->velocity[1]
                << std::setw(8) << atom->velocity[2];

            file << '\n';  // Newline after each ATOM entry

            // Count up iterator
            ++atomCount;

            // Correct for formatting overflow at 100,000 atoms
            if (atomCount == 100000) {
                atomCount = 1;
            }
        }

        // Box vectors
        file << std::fixed << std::setprecision(5)
            << "  " << box_vectors[0] << "  " << box_vectors[1] << "  " << box_vectors[2] << '\n';

        file.close();
    }

}