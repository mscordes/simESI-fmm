#include "SettingsFileReader.h"
#include "Parameters.h"
#include "StringUtils.h" 

#include <fstream> 
#include <iostream>

SettingsFileReader::SettingsFileReader(const fs::path& fileName_) : fileName(fileName_) {};

bool SettingsFileReader::Read() {
    Parameters& params = Parameters::Get();
    ifstream inputFile;
    inputFile.open(fileName.c_str());
    if (inputFile.fail()) throw runtime_error("Could not find the Settings File " + fileName.string() + ".");

    string line;
    while (getline(inputFile, line))
    {
        cleanString(line);
        if (isOnlyWhitespace(line)) continue;

        string arg;
		if (!parseStringBefore(line, arg)) throw runtime_error("Could not parse line: " + line);

        if (!params.getSettingsPresent().test(OUTPUTDIRNAME) && arg == "o") {
            if (parseStringAfter(line, params.getOutputDirName())) {
                params.getOutputDirName() = params.getOutputDirName();
                params.getSettingsPresent().set(OUTPUTDIRNAME);
            }
        }
        else if (!params.getSettingsPresent().test(PDB) && arg == "pdb") {
            if (parseStringAfter(line, params.getPDB())) {
                params.getPDB() = params.getPDB();
                params.getSettingsPresent().set(PDB);
            }
        }
        else if (!params.getSettingsPresent().test(TIME) && arg == "time") {
            if (extractDoubleAfter(line, params.getTime())) {
                params.getSettingsPresent().set(TIME);
            }
        }
        else if (!params.getSettingsPresent().test(GMXENV) && arg == "gmx_env") {
            if (parseStringAfter(line, params.getGMXEnv())) {
                params.getSettingsPresent().set(GMXENV);
            }
        }
        else if (!params.getSettingsPresent().test(GPU) && arg == "gpu") {
            if (parseBoolAfter(line, params.isGPU())) {
                params.getSettingsPresent().set(GPU);
            }
        }
        else if (!params.getSettingsPresent().test(HPC) && arg == "hpc") {
            if (parseBoolAfter(line, params.isHPC())) {
                params.getSettingsPresent().set(HPC);
			}
        }
        else if (!params.getSettingsPresent().test(FMM) && arg == "fmm") {
            if (parseBoolAfter(line, params.isFMM())) {
                params.getSettingsPresent().set(FMM);
            }
        }
        else if (!params.getSettingsPresent().test(VERBOSE) && arg == "verbose") {
            if (parseBoolAfter(line, params.isVerbose())) {
                params.getSettingsPresent().set(VERBOSE);
            }
        }
        else if (!params.getSettingsPresent().test(SAVE) && arg == "save") {
            if (parseBoolAfter(line, params.isSave())) {
                params.getSettingsPresent().set(SAVE);
            }
        }
        else if (!params.getSettingsPresent().test(CONT) && arg == "cont") {
            if (parseBoolAfter(line, params.isCont())) {
                params.getSettingsPresent().set(CONT);
            }
        }
        else if (params.isCont() && !params.getSettingsPresent().test(DIRCONT) && arg == "dir_cont") {
            string path;
            if (parseStringAfter(line, path)) {
                params.getDirCont() = fs::path("outputFiles") / path;
                params.getSettingsPresent().set(DIRCONT);
            }
        }
        else if (params.isCont() && !params.getSettingsPresent().test(STEPCONT) && arg == "step_cont") {
            if (extractIntAfter(line, params.getStepCont())) {
                params.getSettingsPresent().set(STEPCONT);
            }
        }
        else if (!params.getSettingsPresent().test(MODE) && arg == "mode") {
            string mode;
            if (parseStringAfter(line, mode)) {
                if (mode == "pos") params.getMode() = Polarity::POS;
                else if (mode == "neg") params.getMode() = Polarity::NEG;
                else throw runtime_error("Invalid mode in settings file: " + mode + ". Use 'pos' or 'neg'.");
                params.getSettingsPresent().set(MODE);
            }
        }
        else if (!params.getSettingsPresent().test(AMAC) && arg == "amac") {
            if (extractDoubleAfter(line, params.getAmac())) {
                params.getSettingsPresent().set(AMAC);
            }
        }
        else if (!params.getSettingsPresent().test(MPM) && arg == "mpm") {
            if (parseBoolAfter(line, params.isMPM())) {
                params.getSettingsPresent().set(MPM);
            }
        }
        else if (!params.getSettingsPresent().test(DROPLETSIZE) && arg == "droplet_size") {
            if (extractDoubleAfter(line, params.getDropletSize())) {
                params.getSettingsPresent().set(DROPLETSIZE);
            }
        }
        else if (!params.getSettingsPresent().test(WATERCUTOFF) && arg == "water_cutoff") {
            if (extractIntAfter(line, params.getWaterCutoff())) {
                params.getSettingsPresent().set(WATERCUTOFF);
            }
        }
        else if (!params.getSettingsPresent().test(ITEMP) && arg == "itemp") {
            if (extractDoubleAfter(line, params.getITemp())) {
                params.getSettingsPresent().set(ITEMP);
            }
        }
        else if (!params.getSettingsPresent().test(FTEMP) && arg == "ftemp") {
            if (extractDoubleAfter(line, params.getFTemp())) {
                params.getSettingsPresent().set(FTEMP);
            }
        }
        else if (!params.getSettingsPresent().test(PKAPDB) && arg == "pka_pdb") {
            if (parseStringAfter(line, params.getPKAPDB())) {
                params.getPKAPDB() = params.getPKAPDB();
                params.getSettingsPresent().set(PKAPDB);
            }
        }
        else if (!params.getSettingsPresent().test(ATM) && arg == "atm") {
            if (parseBoolAfter(line, params.isATM())) {
                params.getSettingsPresent().set(ATM);
            }
        }
        else if (!params.getSettingsPresent().test(GASTEMP) && arg == "gas_temp") {
            if (extractDoubleAfter(line, params.getGasTemp())) {
                params.getSettingsPresent().set(GASTEMP);
            }
        }
        else if (!params.getSettingsPresent().test(WATERMOLP) && arg == "water_molp") {
            if (extractDoubleAfter(line, params.getWaterMolp())) {
                params.getSettingsPresent().set(WATERMOLP);
            }
        }
        else if (!params.getSettingsPresent().test(N2MOLP) && arg == "n2_molp") {
            if (extractDoubleAfter(line, params.getN2Molp())) {
                params.getSettingsPresent().set(N2MOLP);
            }
        }
        else if (!params.getSettingsPresent().test(O2MOLP) && arg == "o2_molp") {
            if (extractDoubleAfter(line, params.getO2Molp())) {
                params.getSettingsPresent().set(O2MOLP);
            }
        }
        else if (!params.getSettingsPresent().test(BOXSIZE) && arg == "box_size") {
            if (extractDoubleAfter(line, params.getBoxSize())) {
                params.getSettingsPresent().set(BOXSIZE);
            }
        }
        else {
            if (fileName != "inputFiles/default.txt") {
                cerr << "WARNING: Ignored line \"" << line << "\" when parsing " << fileName;
                cerr << " - either unknown line, or parameter has already been set):\n" << endl;
            }
        }
    }

    inputFile.close();
    return true;
}