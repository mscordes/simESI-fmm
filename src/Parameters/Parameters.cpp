#include "Parameters.h"

#include <unordered_map>
#include <stdexcept>

Parameters::Parameters() {
    settingsPresent.reset();
}

Parameters& Parameters::Get() {
    static Parameters instance; 
    return instance;
}

static void printHelp() {
    cout
        << "simESI-fmm - Simulations of ESI with the Fast Multipole Method \n"
        << "  -h or -help          Show this help message and exit.\n"

        << "\nRun Parameters: \n"
        << "  -f <string>           Settings file name in inputFiles/. Default: \"input.txt\".\n"
        << "  -o <string>           (Optional) Set name of output dir in outputFiles/ rather than use .pdb name.\n"
		<< "  -pdb <string>         PDB file name in coordinateFiles/.\n"
        << "  -time <double>        Simulation time in nanoseconds.\n"
		<< "  -gmx_env <string>     GROMACS environment setup command (e.g. \"source /path/to/gromacs/bin/GMXRC\").\n"
		<< "  -gpu <yes|no>         Use GPU acceleration if 'yes', else 'no'.\n"
		<< "  -hpc <yes|no>         Use HPC settings if 'yes', else 'no'.\n"
		<< "  -fmm <yes|no>         Use FMM for long-range electrostatics if 'yes', else 'no'.\n"
		<< "  -verbose <yes|no>     Enable verbose output if 'yes', else 'no'.\n"
		<< "  -save <yes|no>        Save all gmx output files if 'yes', else 'no'.\n"

		<< "\nContinuation Parameters (Optional):\n"
		<< "  -cont <yes|no>        Continue from previous run if 'yes', else 'no'.\n"
		<< "  -dir_cont <string>    Directory to continue from if -cont is 'yes'. Required if continuing.\n"
		<< "  -step_cont <int>      Step to continue from if -cont is 'yes'. Required if continuing.\n"

		<< "\nDroplet Parameters:\n"
		<< "  -mode <pos|neg>       Ionization mode: 'pos' for positively charged droplet or 'neg' for negative.\n"
		<< "  -amac <double>        Droplet ammonium acetate concentration in M.\n"
		<< "  -mpm  <yes|no>        Enable mobile protons for gas phase ion if 'yes', else 'no'.\n"\
		<< "  -droplet_size <double>Initial droplet diameter in nm.\n"
		<< "  -water_cutoff <int>   (Optional) Cutoff for number of water molecules when droplet temperature is ramped from -itemp to -ftemp.\n"
		<< "  -itemp <double>       Initial temperature of the droplet in Kelvin.\n"
		<< "  -ftemp <double>       Final temperature of the droplet in Kelvin.\n"
		<< "  -pka_pdb <string>     (Optional) User different PDB file for pKa calculations.\n"
 
		<< "\nEnvironmental Parameters:\n"
		<< "  -atm <bool>           Simulate ambient atmosphere if 'yes' else 'no' for vacuum.\n"
		<< "  -gas_temp <double>    Gas temperature in Kelvin.\n"
        << "  -water_molp <double>  Water mol% in atmosphere (out of 100). NOTE: -water_molp, -n2_molp, and -o2_molp must add up to 100 +/- 1.\n"
        << "  -n2_molp <double>     Nitrogen mol% in atmosphere (out of 100).\n"
        << "  -o2_molp <double>     Oxygen mol% in atmosphere (out of 100).\n"
        << "  -box_size <double>    Simulation box size in nm.\n"
		<< endl;
}

static bool pathExists(const fs::path& path) {
    if (!filesystem::exists(path)) return false;
    return true;
}

void Parameters::setFromCommandLine(int argc, char** argv) {
    if (argc > 1 && (string(argv[1]) == "-h" || string(argv[1]) == "-help")) {
		printHelp();
		throw runtime_error("Help requested");
    }

    unordered_map<string, string> args;
    for (int i = 1; i < argc; ++i) {
        string key = argv[i];

        if (key[0] != '-') throw runtime_error("Expected flag starting with '-', got: " + key);
        if (i + 1 >= argc) throw runtime_error("Missing value for " + key);

        string value = argv[++i];
        args[key] = value;
    }

    for (const auto& [key, value] : args) {
        if (key == "-f") {
			settingsFile = fs::path("inputFiles") / value;
            settingsPresent.set(SETTINGSFILE);
        }
        else if (key == "-o") {
            outputDirName =  value;
            settingsPresent.set(OUTPUTDIRNAME);
        }
        else if (key == "-pdb") {
            pdb = value;
            settingsPresent.set(PDB);
        }
        else if (key == "-time") {
			this->time = stod(value);
            settingsPresent.set(TIME);
        }
		else if (key == "-gmx_env") {
            gmxEnv = value;
            settingsPresent.set(GMXENV);
        }
        else if (key == "-gpu") {
            if (value == "yes") gpu = true;
            else if (value == "no") gpu = false;
			else throw runtime_error("Invalid value for -gpu: " + value + ". Use 'yes' or 'no'.");
            settingsPresent.set(GPU);
        }
        else if (key == "-hpc") {
            if (value == "yes") hpc = true;
            else if (value == "no") hpc = false;
            else throw runtime_error("Invalid value for -hpc: " + value + ". Use 'yes' or 'no'.");
            settingsPresent.set(HPC);
        }
        else if (key == "-fmm") {
            if (value == "yes") fmm = true;
            else if (value == "no") fmm = false;
            else throw runtime_error("Invalid value for -fmm: " + value + ". Use 'yes' or 'no'.");
            settingsPresent.set(FMM);
        }
        else if (key == "-verbose") {
            if (value == "yes") verbose = true;
            else if (value == "no") verbose = false;
            else throw runtime_error("Invalid value for -verbose: " + value + ". Use 'yes' or 'no'.");
            settingsPresent.set(VERBOSE);
        }
        else if (key == "-save") {
            if (value == "yes") save = true;
            else if (value == "no") save = false;
            else throw runtime_error("Invalid value for -save: " + value + ". Use 'yes' or 'no'.");
            settingsPresent.set(SAVE);
        }
        else if (key == "-cont") {
            if (value == "yes") cont = true;
            else if (value == "no") cont = false;
            else throw runtime_error("Invalid value for -cont: " + value + ". Use 'yes' or 'no'.");
            settingsPresent.set(CONT);
        }
        else if (key == "-dir_cont") {
            dirCont = fs::path("outputFiles") / value;
            settingsPresent.set(DIRCONT);
        }
        else if (key == "-step_cont") {
            stepCont = stoi(value);
            settingsPresent.set(STEPCONT);
        }
        else if (key == "-mode") {
            if (value == "pos") mode = Polarity::POS;
            else if (value == "neg") mode = Polarity::NEG;
            else throw runtime_error("Invalid value for -mode: " + value + ". Use 'pos' or 'neg'.");
            settingsPresent.set(MODE);
        }
        else if (key == "-amac") {
            amac = stod(value);
            settingsPresent.set(AMAC);
        }
        else if (key == "-mpm") {
            if (value == "yes") mpm = true;
            else if (value == "no") mpm = false;
            else throw runtime_error("Invalid value for -mpm: " + value + ". Use 'yes' or 'no'.");
            settingsPresent.set(MPM);
        }
        else if (key == "-droplet_size") {
            dropletSize = stod(value);
            settingsPresent.set(DROPLETSIZE);
        }
        else if (key == "-water_cutoff") {
            waterCutoff = stoi(value);
            settingsPresent.set(WATERCUTOFF);
        }
        else if (key == "-itemp") {
            iTemp = stod(value);
            settingsPresent.set(ITEMP);
        }
        else if (key == "-ftemp") {
            fTemp = stod(value);
            settingsPresent.set(FTEMP);
		}
        else if (key == "-pka_pdb") {
            pkaPdb = value;
        }
        else if (key == "-atm") {
            if (value == "yes") atm = true;
            else if (value == "no") atm = false;
            else throw runtime_error("Invalid value for -atm: " + value + ". Use 'yes' or 'no'.");
			settingsPresent.set(ATM);
        }
        else if (key == "-gas_temp") {
            gasTemp = stod(value);
            settingsPresent.set(GASTEMP);
        }
        else if (key == "-water_molp") {
            waterMolp = stod(value);
            settingsPresent.set(WATERMOLP);
        }
        else if (key == "-n2_molp") {
            n2Molp = stod(value);
            settingsPresent.set(N2MOLP);
        }
        else if (key == "-o2_molp") {
            o2Molp = stod(value);
            settingsPresent.set(O2MOLP);
        }
        else if (key == "-box_size") {
            boxSize = stod(value);
            settingsPresent.set(BOXSIZE);
		}
        else {
            throw runtime_error("Invalid argument: " + key);
        }
    }

	// Use default settings file if not provided
    if (!settingsPresent.test(SETTINGSFILE)) { 
        settingsFile = fs::path("inputFiles") / "default.txt";
		settingsPresent.set(SETTINGSFILE);
    }
}

void Parameters::areValid() {
	if (!settingsPresent.test(SETTINGSFILE)) throw runtime_error("Missing -f for settings file.");
	if (!settingsPresent.test(PDB)) throw runtime_error("Missing -pdb for PDB file.");
	if (!settingsPresent.test(TIME)) throw runtime_error("Missing -time for simulation time.");
	if (!settingsPresent.test(GMXENV)) throw runtime_error("Missing -gmx_env for GROMACS environment.");
	if (!settingsPresent.test(GPU)) throw runtime_error("Missing -gpu for GPU acceleration.");
	if (!settingsPresent.test(HPC)) throw runtime_error("Missing -hpc for HPC settings.");
	if (!settingsPresent.test(FMM)) throw runtime_error("Missing -fmm for FMM usage.");
	if (!settingsPresent.test(VERBOSE)) throw runtime_error("Missing -verbose for verbose output.");
	if (!settingsPresent.test(SAVE)) throw runtime_error("Missing -save for saving output files.");
	if (!settingsPresent.test(CONT)) throw runtime_error("Missing -cont for continuation.");
	if (!settingsPresent.test(MODE)) throw runtime_error("Missing -mode for ionization mode.");
	if (!settingsPresent.test(AMAC)) throw runtime_error("Missing -amac for ammonium acetate concentration.");
	if (!settingsPresent.test(MPM)) throw runtime_error("Missing -mpm for mobile protons.");
	if (!settingsPresent.test(DROPLETSIZE)) throw runtime_error("Missing -droplet_size for droplet size.");
	if (!settingsPresent.test(ITEMP)) throw runtime_error("Missing -itemp for initial temperature.");
	if (!settingsPresent.test(FTEMP)) throw runtime_error("Missing -ftemp for final temperature.");
	if (!settingsPresent.test(ATM)) throw runtime_error("Missing -atm for atmosphere.");
	if (!settingsPresent.test(GASTEMP)) throw runtime_error("Missing -gas_temp for gas temperature.");
	if (!settingsPresent.test(WATERMOLP)) throw runtime_error("Missing -water_molp for water mol%.");
	if (!settingsPresent.test(N2MOLP)) throw runtime_error("Missing -n2_molp for nitrogen mol%.");
	if (!settingsPresent.test(O2MOLP)) throw runtime_error("Missing -o2_molp for oxygen mol%.");
	if (!settingsPresent.test(BOXSIZE)) throw runtime_error("Missing -box_size for box size.");

    if (!pathExists(settingsFile)) throw runtime_error("Could not find " + settingsFile.string() + " for -f.");
    if (!pathExists("coordinateFiles/" + pdb)) throw runtime_error("Could not find " + pdb + " for -pdb.");

    if (settingsPresent.test(PKAPDB) && !pathExists("coordinateFiles/" + pkaPdb))
        throw runtime_error("Could not find " + pkaPdb + " for -pka_pdb.");

    if (cont) {
        if (!settingsPresent.test(DIRCONT)) throw runtime_error("Missing -dir_cont for continuation.");
        else if (!pathExists(dirCont)) throw runtime_error("Could not find " + dirCont.string() + " for -dir_cont.");

        if (!settingsPresent.test(STEPCONT)) throw runtime_error("Missing -step_cont for continuation.");
        else if (stepCont < 0) throw runtime_error("Invalid -step_cont: " + to_string(stepCont) + ". Must be >= 0.");
    }

    if (time < 0) throw runtime_error("Invalid simulation time. Must be greater than 0.");
    if (amac < 0) throw runtime_error("Invalid ammonium acetate concentration. Must be greater than or equal to 0.");
    if (dropletSize < 0) throw runtime_error("Invalid droplet size. Must be greater than 0.");
    if (iTemp < 0)  throw runtime_error("Invalid initial temperature. Must be greater than 0.");
    if (fTemp < 0) throw runtime_error("Invalid final temperature. Must be greater than 0.");
    if (gasTemp < 0) throw runtime_error("Invalid gas temperature. Must be greater than 0.");
    if (boxSize < 0) throw runtime_error("Invalid box size. Must be greater than or equal to 0.");
    if (waterMolp < 0) throw runtime_error("Invalid water molp. Must be greater than 0.");
    if (n2Molp < 0) throw runtime_error("Invalid N2 molp. Must be greater than 0.");
    if (o2Molp < 0) throw runtime_error("Invalid O2 molp. Must be greater than 0.");
	if (waterCutoff < 0) throw runtime_error("Invalid water cutoff. Must be greater than or equal to 0.");

    double tot_molp = waterMolp + n2Molp + o2Molp;
    if (tot_molp < 99 || tot_molp > 101)
        throw runtime_error("Invalid gas composition. Total mol% must be 100 +/- 1. Current total: " + to_string(tot_molp));
}

void Parameters::print(ostream& os) const {
    os
        << "# Parameters\n"
        << "  o                 = " << outputDirName << "\n"
		<< "  pdb               = " << pdb << "\n"
		<< "  time              = " << time << "\n"
		<< "  gmx_env           = " << gmxEnv << "\n"
		<< "  gpu               = " << (gpu ? "yes" : "no") << "\n"
		<< "  hpc               = " << (hpc ? "yes" : "no") << "\n"
		<< "  fmm               = " << (fmm ? "yes" : "no") << "\n"
		<< "  verbose           = " << (verbose ? "yes" : "no") << "\n"
		<< "  save              = " << (save ? "yes" : "no") << "\n"
		<< "  cont              = " << (cont ? "yes" : "no") << "\n"
		<< "  dir_cont          = " << (cont ? dirCont : "N/A") << "\n"
		<< "  step_cont         = " << (cont ? to_string(stepCont) : "N/A") << "\n"
		<< "  mode              = " << (mode == Polarity::POS ? "pos" : "neg") << "\n"
		<< "  amac              = " << amac << "\n"
		<< "  mpm               = " << (mpm ? "yes" : "no") << "\n"
		<< "  droplet_size      = " << dropletSize << "\n"
		<< "  water_cutoff      = " << waterCutoff << "\n"
		<< "  itemp             = " << iTemp << "\n"
		<< "  ftemp             = " << fTemp << "\n"
		<< "  pka_pdb           = " << (settingsPresent.test(PKAPDB) ? pkaPdb : "N/A") << "\n"
		<< "  atm               = " << (atm ? "yes" : "no") << "\n"
		<< "  gas_temp          = " << gasTemp << "\n"
		<< "  water_molp        = " << waterMolp << "\n"
		<< "  n2_molp           = " << n2Molp << "\n"
		<< "  o2_molp           = " << o2Molp << "\n"
        << "  box_size          = " << boxSize << "\n"
		<< endl;
}

bitset<SET_COUNT>& Parameters::getSettingsPresent() {
	return settingsPresent;
}

fs::path& Parameters::getSettingsFile() {
    return settingsFile;
}

string& Parameters::getOutputDirName() {
    return outputDirName;
}

string& Parameters::getPDB() {
    return pdb;
}

double& Parameters::getTime() {
    return time;
}

string& Parameters::getGMXEnv() {
    return gmxEnv;
}

bool& Parameters::isGPU() {
    return gpu;
}
bool& Parameters::isHPC() {
    return hpc;
}

bool& Parameters::isFMM() {
    return fmm;
}

bool& Parameters::isVerbose() {
    return verbose;
}

bool& Parameters::isSave() {
    return save;
}

bool& Parameters::isCont() {
    return cont;
}

fs::path& Parameters::getDirCont() {
    return dirCont;
}

int& Parameters::getStepCont() {
    return stepCont;
}

Polarity& Parameters::getMode() {
    return mode;
}

double& Parameters::getAmac() {
    return amac;
}

bool& Parameters::isMPM() {
    return mpm;
}

double& Parameters::getDropletSize() {
    return dropletSize;
}

int& Parameters::getWaterCutoff() {
    return waterCutoff;
}

double& Parameters::getITemp() {
    return iTemp;
}

double& Parameters::getFTemp() {
    return fTemp;
}

string& Parameters::getPKAPDB() {
    return pkaPdb;
}

bool& Parameters::isATM() {
    return atm;
}

double& Parameters::getGasTemp() 
{ 
    return gasTemp; 
}

double& Parameters::getWaterMolp() 
{ 
    return waterMolp; 
}

double& Parameters::getN2Molp() 
{ 
    return n2Molp; 
}

double& Parameters::getO2Molp() 
{ 
    return o2Molp; 
}

double& Parameters::getBoxSize() 
{ 
    return boxSize; 
}

bool Parameters::isHumid() const
{
    return (waterMolp > 0.001);
}

fs::path& Parameters::getOutputDir() {
    return outputDir;
}

fs::path& Parameters::getDataDir() {
    return dataDir;
}

fs::path& Parameters::getSimulationDir() {
    return simulationDir;
}