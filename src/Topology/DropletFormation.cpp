#include "DropletFormation.h"
#include "Atom.h"
#include "Atmosphere.h"
#include "Constants.h"
#include "CoordinateManip.h"
#include "Charge.h"
#include "FileUtils.h"
#include "Parameters.h"
#include "TopologyUtils.h"
#include "Subprocess.h"

#include <cmath>
#include <iostream>
#include <numeric>

// Complete droplet formation, note that droplet is formed as a shell around the protein
void formDroplet(State& state) {
	cout << "Beginning droplet formation." << endl;

	Parameters& params = Parameters::Get();
	const string& gmxEnv = params.getGMXEnv();

	const auto& boxSize = params.getBoxSize();
	Vec3D boxD(boxSize, boxSize, boxSize);

	// Make box as small as possible to ease water seeding
	Vec3D smallBoxD = centerSmallBox(state); 
	Vec3D center = smallBoxD / 2;
	state.writeGRO("pre_solvate.gro", smallBoxD);

	vector<Vec3D> protCoords = getProteinCoords(state);

	// Find number of waters in droplet for calculating solutes to seed
	int numWater;
	{ 
		copyFile("system.top", "temp.top");
		string command =  gmxEnv + " solvate -cp pre_solvate.gro -cs tip4p_2005.gro -p temp.top -o temp.gro";
		subprocess(command);

		State tempState;
		tempState.parseGRO("temp.gro");
		carveDroplet(tempState);
		numWater = tempState.residueSet.at("SOL").size();
	}

	// Get radius of droplet assuming sphereical vol
	double protVol = getProteinMass(state) / (Constants::NA * 1000 * Constants::protDensity); 
	double waterVol = numWater * Constants::waterVol;
	double dropletVol = protVol + waterVol;
	double dropletR = cbrt((3 * dropletVol) / (4 * Constants::pi)); // m

	// Calculate 90% of Rayleigh limit as initial droplet net charge
	double rayleighLim = Constants::rayleighC * pow(dropletR, 3.0 / 2.0);
	int dropletCharge = static_cast<int>(round(0.9 * rayleighLim));

	// Get charge of solute to seed accounting for protein charge
	int protCharge = static_cast<int>(round(getProteinCharge(state)));
	int soluteCharge = dropletCharge - abs(protCharge);

	// Begin finding how many of each residue type to add now that we have number of waters
	unordered_map<string, int> tooSeed;
	double numAmAc = 0.01801 * numWater * params.getAmac();

	// Excess AmAc dependent on ESI mode (for net droplet charge), also effect from pH
	Polarity mode = params.getMode();
	if (mode == Polarity::POS) {
		tooSeed["ATX"] = static_cast<int>(round(numAmAc * 0.9));	
		tooSeed["AHX"] = static_cast<int>(round(numAmAc * 0.1));	
		tooSeed["NXH"] = soluteCharge + tooSeed["ATX"];		
	}
	else if (mode == Polarity::NEG) {
		tooSeed["NXH"] = static_cast<int>(round(numAmAc * 0.9));	
		tooSeed["NXX"] = static_cast<int>(round(numAmAc * 0.1));
		tooSeed["ATX"] = soluteCharge + tooSeed["NXH"];		
	}
	else throw runtime_error("Invalid ESI mode.");

	// For no droplet
	if (abs(params.getDropletSize()) < numeric_limits<double>::epsilon())
		for (auto& pair : tooSeed)
			tooSeed[pair.first] = 0;

	vector<Vec3D> coords = protCoords;
	insertMolecIntoDroplet("ATX", tooSeed, state, protCoords, coords, smallBoxD);
	insertMolecIntoDroplet("AHX", tooSeed, state, protCoords, coords, smallBoxD);
	insertMolecIntoDroplet("NXH", tooSeed, state, protCoords, coords, smallBoxD);
	insertMolecIntoDroplet("NXX", tooSeed, state, protCoords, coords, smallBoxD);

	{ // Resolvate + atmosphere
		state.writeGRO("solute.gro", smallBoxD);
		string command = gmxEnv +  " solvate -cp solute.gro -cs tip4p_2005.gro -p temp.top -o solvated.gro";
		subprocess(command);

		state.parseGRO("solvated.gro");
		carveDroplet(state); 

		centerState(state, boxD);
		if (params.isATM()) seedAtmosphere(state);
	}

	createTOP("droplet.top", state);
	state.writeGRO("droplet.gro", boxD);
	cout << "Droplet formation completed." << endl;
}