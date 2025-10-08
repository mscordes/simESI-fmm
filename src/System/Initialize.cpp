#include "System.h"
#include "Atom.h"
#include "CoordinateManip.h"
#include "DropletFormation.h"
#include "FileUtils.h"
#include "GmxUtils.h"
#include "ParsePKA.h"
#include "Protein.h"
#include "Residue.h"
#include "Subprocess.h"
#include "TopologyUtils.h"

#include <iostream>
#include <fstream>

// Process inputted .pdb, form droplet, then equilibriate 
void System::initialize() {
    const string& pdb = params.getPDB();
    const string& gmxEnv = params.getGMXEnv();

    // Clean .pdb of any non-protein residues and center system
    if (!fs::exists(pdb)) throw runtime_error("Could not find inputted .pdb file.");
    state.parseProteinPDB(pdb);
    centerState(state, params.getBoxSize());
    state.writeProteinPDB("init.pdb");

    { // Test validity of inputted .pdb and convert to a properly formatted .gro
        vector<string> inputs;
        for (const auto& [name, chain] : state.proteinMap) {
            string nterResName = chain->residues.front()->name;
            inputs.push_back((nterResName == "PRO" || nterResName == "MET") ? "1" : "0"); // Nterm code special for MET/PRO
            inputs.push_back("0");
        }

        // pdb2gmx
        string command = gmxEnv + " pdb2gmx -f init.pdb -o init.gro -p temp.top -ff charmm36 -water none -ignh -ter";
        subprocess(command, inputs);

        if (!filesystem::exists("temp.top"))
            throw runtime_error("Could not form molecular topology. Ensure 'gmx' is callable.");

        string dummy;
        if (!(ifstream("temp.top") >> dummy))
            throw runtime_error("Could not form molecular topology. Input .pdb likely compromised.");
    }

    // Update State with now properly formatted .gro
    state.parseGRO("init.gro");

    // Make "system.gro" and "system.top" files with proper protonation states and update state
    parsePKA(state);
    double ph = (Parameters::Get().getMode() == Polarity::POS) ? 5.5 : 8.5;

    call_pdb2gmx("system.top", setPdb2GmxInputs(state, ph), state, true);
    state.parseGRO("system.gro");

    // Form droplet and atmosphere
    formDroplet(state);

    { // Energy minimization
        string command = gmxEnv + " grompp -f em.mdp -c droplet.gro -p droplet.top -o em.tpr -maxwarn 100";
        subprocess(command);
        call_mdrun("em");
    }

    { // NVT
        modifyMDPgrps("nvt", state, Parameters::Get().getITemp());
        string command = gmxEnv + " grompp -f nvt.mdp -c em.gro -p droplet.top -n nvt.ndx -o nvt.tpr -maxwarn 100";
        subprocess(command);
        call_mdrun("nvt");
    }

    if (fs::exists("nvt.gro")) cout << "Equilibriation completed successfully.\n";
    else throw runtime_error("Equilibriation failed.");
    state.parseGRO("nvt.gro");

    copyFile("nvt.gro", "0.gro");
    copyFile("droplet.top", "0.top");
}