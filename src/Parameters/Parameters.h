#pragma once
#include <bitset>
#include <filesystem>
#include <iostream>
#include <string>

using namespace std;
namespace fs = filesystem;

enum class Polarity {
	POS,
	NEG
};

enum SettingsFileContents {
	SETTINGSFILE,
	OUTPUTDIRNAME,
	PDB,
	TIME,
	GMXENV,
	GPU,
	HPC,
	FMM,
	VERBOSE,
	SAVE,
	CONT, 
	DIRCONT, 
	STEPCONT,
	MODE, 
	AMAC,
	MPM,
	DROPLETSIZE,
	WATERCUTOFF,
	ITEMP,
	FTEMP,
	PKAPDB,
	ATM, 
	GASTEMP,
	WATERMOLP,
	N2MOLP,
	O2MOLP,
	BOXSIZE,

	SET_COUNT   // sentinel, always keep last
};

class Parameters {
	private:
		static Parameters*	instance;
		bitset<SET_COUNT>	settingsPresent{};

		// Run Parameters
		fs::path			settingsFile;
		string				outputDirName;
		string				pdb;
		double				time;
		string				gmxEnv;
		bool				gpu;
		bool				hpc;
		bool				fmm;
		bool				verbose;
		bool				save;

		// Continuation Parameters
		bool				cont;
		fs::path			dirCont;
		int 				stepCont;

		// Droplet Parameters
		Polarity			mode;
		double 				amac;
		bool 				mpm;
		double              dropletSize;
		int					waterCutoff;
		double				iTemp;
		double				fTemp;
		string				pkaPdb;

		// Environmental Parameters
		bool 				atm;
		double				gasTemp;
		double              waterMolp;
		double              n2Molp;
		double              o2Molp;
		double              boxSize;

		// Output directories
		fs::path			outputDir;
		fs::path 			dataDir;	
		fs::path 			simulationDir;

	public:
		static Parameters & Get();
							Parameters();   
		void				setFromCommandLine(int argc, char** argv);
		void				areValid();
		void			    print(ostream& os) const;

		bitset<SET_COUNT> & getSettingsPresent();
		fs::path		  & getSettingsFile();
		string			  & getOutputDirName();
		string			  & getPDB();
		double			  & getTime();
		string			  & getGMXEnv();
		bool			  & isGPU();
		bool			  & isHPC();
		bool			  & isFMM();
		bool			  & isVerbose();
		bool			  & isSave();
		bool			  & isCont();
		fs::path		  & getDirCont();
		int				  & getStepCont();
		Polarity	      & getMode();
		double			  & getAmac();
		bool			  & isMPM();
		double			  & getDropletSize();
		int				  & getWaterCutoff();
		double			  & getITemp();
		double			  & getFTemp();
		string			  & getPKAPDB();
		bool			  & isATM();
		double			  & getGasTemp();
		double			  & getWaterMolp();
		double			  & getN2Molp();
		double			  & getO2Molp();
		double			  & getBoxSize();
		bool			    isHumid() const;
		fs::path		  & getOutputDir();
		fs::path		  & getDataDir();
		fs::path		  & getSimulationDir();
};