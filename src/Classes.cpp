#include <iostream>
#include <string>
#include <array> 
#include <vector> 
#include <map>
#include <unordered_map>
#include <math.h>
#include <memory>

namespace Core {

    // User arguments
    struct Config {
        std::string pdb;
        std::string esi_mode = "pos";
		std::string atm = "yes";
        float amace_conc = 0.25f;
        float water_vapor = 0.00f;
        float ace_vapor = 0.00f;
        float ach_vapor = 0.00f;
        float nh4_vapor = 0.00f;
        float nh3_vapor = 0.00f;
        float droplet_size = 1.1f;
        float box_size = 75.0f;
        float time = 25.0f;
        float init_temp = 370.0f;
        float final_temp = 450.0f;
        float gas_temp = 300.0f;
		int water_cutoff = -1;
        std::string pka_pdb = "no";
        std::string dir_cont = "no";
        int step_cont = -1;
		std::string gmx_env = "gmx";
		std::string gpu = "yes";
        std::string hpc = "no";
        std::string fmm = "no";
    };

    struct Residue; // Forward declaration of Residue

    struct Atom {
    public:
        // Attributes
		std::weak_ptr<Residue> parent;      // Pointer to parent Residue object
        std::string res_id;                 // Residue number + residue name
        int res_num;                        // Residue number
        std::string res_name;               // Residue name
        std::string atom_name;              // Atom name
        int atom_num;                       // Atom number
        std::array<float, 3> coord;         // 3D coordinate of atom
        std::array<float, 3> velocity;      // 3D atom velocity
        std::string element;                // Atom element
        std::string chain;                  // Chain

        // Default constructor (required for containers)
        Atom() = default;

        Atom(std::shared_ptr<Residue> parent, std::string res_id, int res_num, 
            std::string res_name, std::string atom_name, int atom_num,
            std::array<float, 3> coord, std::array<float, 3> velocity,
            std::string element, std::string chain)
            : res_id(std::move(res_id)), res_num(res_num), res_name(std::move(res_name)),
            atom_name(std::move(atom_name)), atom_num(atom_num), coord(std::move(coord)),
            velocity(std::move(velocity)), element(std::move(element)), chain(std::move(chain)) {}

        // String representation 
        std::string toString() const {
            return "Atom: " + atom_name +
                " Residue: " + res_id + " Chain: " 
                + chain + " Coordinates: [" +
                std::to_string(coord[0]) + ", " +
                std::to_string(coord[1]) + ", " +
                std::to_string(coord[2]) + "]";
        }
    };

    struct Residue {
    public:
        // Attributes
        std::string res_id;                         // Residue number + residue name
        int res_num;                                // Residue number
        std::string res_name;                       // Residue name
        std::string chain;                          // Chain
        std::vector<std::weak_ptr<Atom>> atoms;   // Pointers to Atom objects that compose the residue

        // Default constructor (required for containers)
        Residue() : res_id(""), res_num(0), res_name(""), chain(""), atoms() {}

        Residue(std::string res_id, int res_num, std::string res_name, std::string chain)
            : res_id(std::move(res_id)), res_num(res_num), res_name(std::move(res_name)), chain(std::move(chain)) {}

        // Method to add an Atom to the Residue
        void addAtom(const std::weak_ptr<Atom>& atom) {
            atoms.push_back(atom); // Use a shared pointer to the original atom
        }

        // String Representation
        std::string toString() const {
            std::string result = "Residue: " + res_id + " Chain: " + chain + "\nAtoms:\n";
            for (const auto& atom : atoms) {
                result += "  " + atom.lock()->toString() + "\n";
            }
            return result;
        }
    };

    struct Protein {
    public:
        // Attributes
        std::string chain;                              // Chain
        std::vector<std::weak_ptr<Residue>> residues; // Pointers to Residue objects that compose the protein

        // Default constructor (required for containers)
        Protein() = default;

        // Constructor
        Protein(std::string chain) : chain(std::move(chain)) {}

        // Method to add a Residue to the Protein
        void addResidue(const std::weak_ptr<Residue>& residue) {
            residues.push_back(residue); 
        }

        // String Representation
        std::string toString() const {
            std::string result = "Chain: " + chain + "\nFrom Residues: ";
            if (!residues.empty()) {
                result += residues[0].lock()->res_id;
                result += " to ";
                result += residues.back().lock()->res_id;
            }
            result += ".";
            return result;
        }
    };

    // Combined object with condensed information of a particular coordinate file
    struct CoordInfo {
    public:
        // Attributes
        std::vector<std::shared_ptr<Atom>> atoms;                               // Atoms
        std::vector<std::shared_ptr<Residue>> residues;                         // Residues
        std::vector<std::shared_ptr<Protein>> proteins;                         // Proteins
        std::map<std::string, std::vector<std::shared_ptr<Atom>>> proteinAtoms; // Protein atoms
        std::vector<std::array<float, 3>> coordinates;                          // Coordinates
        std::vector<std::array<float, 3>> proteinCoords;                        // Protein atom coordinates
		std::vector<std::array<float, 3>> waterOCoords;                         // Water O coordinates 
		std::vector<std::array<float, 3>> waterHCoords;                         // Water H coordinates
        std::array<float, 3> box_vectors;                                       // Simulation box vectors
        std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>> residueMap; // Residue Map
        std::unordered_map<std::string, int> numResidues;                       // Number of each residue type

        CoordInfo(
            std::vector<std::shared_ptr<Atom>> atoms = {},
            std::vector<std::shared_ptr<Residue>> residues = {},
            std::vector<std::shared_ptr<Protein>> proteins = {},
            std::map<std::string, std::vector<std::shared_ptr<Atom>>> proteinAtoms = {},
            std::vector<std::array<float, 3>> coordinates = {},
            std::vector<std::array<float, 3>> proteinCoords = {},
            std::vector<std::array<float, 3>> waterOCoords = {},
			std::vector<std::array<float, 3>> waterHCoords = {},
            std::array<float, 3> box_vectors = { 0.0f, 0.0f, 0.0f },
            std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>> residueMap = {},
            std::unordered_map<std::string, int> numResidues = {}
        ) : atoms(std::move(atoms)),
            residues(std::move(residues)),
            proteins(std::move(proteins)),
            proteinAtoms(std::move(proteinAtoms)),
            coordinates(std::move(coordinates)),
            proteinCoords(std::move(proteinCoords)),
			waterOCoords(std::move(waterOCoords)),
			waterHCoords(std::move(waterHCoords)),
            box_vectors(box_vectors),
            residueMap(std::move(residueMap)),
            numResidues(std::move(numResidues))
        {}
    };

    struct Cluster {
        int clusterID;
        int numWaters;
        int numH3O;
		int numOH;
        float pH;

        // String Representation
        std::string toString() const {
            std::string result = 
                "Cluster Id: " + std::to_string(clusterID) + 
                "\nWaters  : " + std::to_string(numWaters) + 
                "\nH3O+    : " + std::to_string(numH3O) +
                "\nOH-     : " + std::to_string(numOH) +
				"\npH      : " + std::to_string(pH);
            return result;
        }
    };

	// Information required to compute exchange for a particular titratable atom (either proton donor or acceptor)
    struct TitratableAtom {
    public:
        // Attributes
        std::shared_ptr<Atom> atom;  // Atom object
        std::array<float, 3> coord;  // Atom coordinates
        int nearWaters;              // Number of coordinated waters
        int clusterID;               // Cluster ID of cluster solvating atom
        float pka;                   // pKa of the titratable atom  
        float gpb;                   // Gas phase basicity (kJ/mol)         

        TitratableAtom(
            std::shared_ptr<Atom> atom, std::array<float, 3> coord, int nearWaters, int clusterID, float pka, float gpb
        ) : atom(atom), coord(coord), nearWaters(nearWaters), clusterID(clusterID), pka(pka), gpb(gpb)
        {}

        // String representation 
        std::string toString() const {
            return atom->toString() + " Near Waters: " + std::to_string(nearWaters) +
                " Cluster ID: " + std::to_string(clusterID) + " pKa: " + std::to_string(pka) +
                " GPB: " + std::to_string(gpb);
        }
    };

	// Combine all titratable atoms into a single object
    struct TitratableSites {
    public:
        // Attributes
        std::vector<TitratableAtom> donors;     // Proton donors
		std::vector<TitratableAtom> acceptors;  // Proton acceptors

		TitratableSites(
			std::vector<TitratableAtom> donors = {},
			std::vector<TitratableAtom> acceptors = {}
		) : donors(std::move(donors)),
			acceptors(std::move(acceptors)) 
        {}

        // String Representation
        std::string toString() const {
            std::string result = "Donors\n";
            for (const auto& atom : donors) {
                result += "  " + atom.toString() + "\n";
            }
			result += "\n\nAcceptors\n";
			for (const auto& atom : acceptors) {
				result += "  " + atom.toString() + "\n";
			}
            return result;
        }
    };

    // All information related to a particular exchange
    struct Exchange {
    public:
        // Attributes
        std::shared_ptr<Atom> donor;     // Hydrogen donor
        std::shared_ptr<Atom> acceptor;  // Accepting site
        float energy;   // Energy of the exchange
        int step;	    // Step at which exchange occurred
        int hop;        // Hop at which exchange occurred
        int clusterID;  // Cluster ID of cluster solvating pair
		float ph;       // pH of solvating cluster
        int nearWaters; // Number waters in solvating cluster      

        Exchange(
            std::shared_ptr<Atom> donor, std::shared_ptr<Atom> acceptor, float energy, int step, 
            int hop, int clusterID, float ph, int nearWaters
		) : donor(donor), acceptor(acceptor), energy(energy), step(step), hop(hop), clusterID(clusterID), 
            ph(ph), nearWaters(nearWaters)
        {}

        // String representation 
        std::string toString() const {
            return "Donor: " + donor->toString() + "\nAcceptor: " + acceptor->toString() + "\nEnergy: " + std::to_string(energy) +
                "\nStep: " + std::to_string(step) + "\nHop: " + std::to_string(hop) + "\nCluster ID: " + std::to_string(clusterID) +
                "\npH: " + std::to_string(ph) + "\nNear Waters: " + std::to_string(nearWaters) + "\n";
        }
    };
}