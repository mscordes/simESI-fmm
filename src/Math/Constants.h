#pragma once
#include <cmath>
#include <limits>
#include <numbers>
#include <unordered_map>
#include <string>
#include <vector>
#include <unordered_set>

using namespace std;

namespace Constants {

    // Physical constants
    constexpr double inf = numeric_limits<double>::infinity();
    constexpr double pi = numbers::pi;
	constexpr double e = 1.602176634e-19;           // Elementary charge (C)
    constexpr double eo = 8.8541878188e-12;         // Vacuum permittivity (F/m)
    constexpr double kB = 1.380649e-23;             // Boltzmann constant (J/K)
    constexpr double R = 8.31446261815324;          // Gas constant (J/mol*K)
    constexpr double NA = 6.02214076e23;            // Avogadros number
    constexpr double ke = 138.932;                  // Coulomb constant (kJ * nm / mol * e^2)
    constexpr double standardAtm = 101325;          // Standard atmosphere (Pa)
    constexpr double waterVol = 2.989040186e-29;    // Volume of single water molecule (m^3)
    constexpr double protDensity = 1220;            // Protein density (kg/m^3) 
    constexpr double surfaceT = 0.0728;             // Water surface tension (N/m)
    const double rayleighC = ((8 * pi) / e) * sqrt(eo * surfaceT); // Reduced constant for Rayleigh Limit calc's

    // M for virtual site
    inline const unordered_map<char, double> atomMasses = {
        {'H', 1.008}, {'C', 12.011}, {'N', 14.007}, {'O', 15.999}, {'S', 32.06}, {'M', 0.000}
    };

    inline const unordered_map<string, double> defaultPkaVals = {
        { "LYS", 10.54 }, { "ARG", 12.48 }, { "ASP", 3.90 }, { "GLU", 4.07 }, { "HIS", 6.04 },
        { "N+", 9.00 }, { "C-", 2.00 }, {"HHO", 0.0}, {"OHX", 14.0}, {"ATX", 4.76}, {"AHX", 4.76}, 
        {"NXX", 9.25}, {"NXH", 9.25}
    };

    // Gas phase basicities (kJ/mol)
    inline const unordered_map<string, double> GPBs = {
        {"HIS", 935.54}, {"LYS", 884.08}, {"ARG", 983.24}, {"GLU", 1424.23}, {"ASP", 1428.42}, 
        {"NTERM", 850.90}, {"CTERM", 1400.38}, {"OHX", 1605.40}, {"HHO", 659.82}, {"ATX", 1428.42}, 
        {"AHX", 1428.42}, {"NXX", 818.81}, {"NXH", 818.81}
    };

	inline const string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    inline const unordered_map<string, string> templateGroNames = {
        {"SOL", "sol.gro"}, {"HHO", "h3o.gro"}, {"OHX", "oh.gro"}, {"ATX", "ace.gro"}, {"AHX", "aceh.gro"}, 
        {"NXX", "nh3.gro"}, {"NXH", "nh4.gro"}, {"NNN", "n2.gro"}, {"OOO", "o2.gro"}
    };

    inline const unordered_set<string> aminoAcids = {
        "ALA","ARG","ASN","ASP","CYS",
        "GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO",
        "SER","THR","TRP","TYR","VAL"
    };

    // All titratable residues excluding water
    inline const unordered_set<string> titAminoAcids = {
        "LYS", "ARG", "ASP", "GLU", "HIS"
    };

    // Order of residues in .top/.gro files
    inline const vector<string> topOrder = { "SOL", "HHO", "OHX", "ATX", "AHX", "NXX", "NXH", "NNN", "OOO" };

    // Names of atoms that can be last in a chain
    inline const unordered_set<string> finalAtoms = { "OT2", "HT2", "OXT", };

    // Order that pdb2gmx takes residues (exclude termini here)
    inline const vector<string> pdb2gmxOrder = { "LYS", "ARG", "ASP", "GLU", "HIS" };

    // All titratable residues excluding water
    inline const unordered_set<string> titResidues = {
        "LYS", "ARG", "ASP", "GLU", "HIS", "HHO", "OHX", "ATX", "AHX", "NXX", "NXH"
    };

    // All titratable residues excluding water and protein
    inline const unordered_set<string> nonProtTitResidues = {
        "HHO", "OHX", "ATX", "AHX", "NXX", "NXH"
    };
    
    // Titrable hydrogens of each titrable amino acid excluding water (CHARMM36)
    inline const unordered_map<string, unordered_set<string>> hydrogensMap = {
        { "LYS", {"HZ1", "HZ2", "HZ3"}}, { "ARG", {"HH11", "HH12", "HH21", "HH22"}},
        { "ASP", {"HD2"}}, { "GLU", {"HE2"}}, { "HIS", {"HD1", "HE2"}},
        {"NTERM", {"H1", "H2", "H3"}},  {"CTERM", {"HT2"}}, 
        { "HHO", {"HW1", "HW2", "HW3"}}, {"AHX", {"HO1"}}, {"NXH", {"HZ1", "HZ2", "HZ3", "HZ4"}}
    };

    // Hydrogen acceptor names of each titrable amino acid excluding water (CHARMM36)
    inline const unordered_map<string, unordered_set<string>> acceptorsMap = {
        { "LYS", {"NZ"} }, { "ARG", {"NH1"} }, { "ASP", {"OD1", "OD2"}}, {"GLU", {"OE1", "OE2"}}, {"HIS", {"NE2", "ND1"}},
        {"NTERM", {"N"}},  {"CTERM", {"OT1", "OT2"}},
        { "OHX", {"O1"}}, {"ATX", {"O1", "O2"}}, {"NXX", {"N1"}}
    };

    // pdb2gmx inputs for titratable residues (excluding N-termini), first = prot, second = deprot
    // Assume HISD here for now, HISE/HISD set seperately as this is meant for initial prot states
    inline const unordered_map<string, pair<string, string>> pdb2gmxCodes = {
        { "LYS", {"1", "0"} }, { "ARG", {"1", "0"} }, { "ASP", {"1", "0"} },
        { "GLU", {"1", "0"} }, { "HIS", {"2", "0"} }, { "CTERM", {"1", "0"} }
    };

    // Partial charges of non-protein residues (CHARMM36)
    inline const unordered_map<string, vector<double>> partialsMap = {
        { "SOL", {    0., 0.5564, 0.5564, -1.1128 }},
        { "HHO", {    0., 0.416 , 0.416 , -0.248 , 0.416 }},
        { "OHX", { -1.32, 0.32 }},
        { "NNN", {    0., 0. }},
        { "OOO", {    0., 0. }},
        { "ATX", { -0.37, 0.62  , 0.09  , 0.09  , 0.09  , -0.76 , -0.76 }},
        { "AHX", { -0.3 , 0.75  , 0.09  , 0.09  , 0.09  , -0.55 , -0.6  , 0.43 }},
        { "NXX", { -1.125, 0.375 , 0.375 , 0.375 }},
        { "NXH", { 0.33 , -0.33 , 0.33  , 0.33  , 0.33 }}
    };

}