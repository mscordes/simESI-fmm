#include "Core.h"		
#include "Parameters.h"	
#include "ParamInfoWriter.h"
#include "FileUtils.h"
#include "SettingsFileReader.h"
#include "System.h"

#include <iostream>
#include <filesystem>

using namespace std;
namespace fs = filesystem;

void run(int argc, char** argv) {

	// Initialize parameters
	Parameters& params = Parameters::Get(); 
	params.setFromCommandLine(argc, argv);
	SettingsFileReader r(params.getSettingsFile());
	r.Read();
	SettingsFileReader rDefault(fs::path("inputFiles") / "default.txt"); 
	rDefault.Read();
	params.areValid();
	params.print(cout);

	// Initialize output directories
	populateOutputDirs(); 
	ParamInfoWriter pWriter(params.getOutputDir() / "parameters.txt");
	pWriter.write();
	fs::current_path(params.getSimulationDir());

	System system; 
	system.run();
}