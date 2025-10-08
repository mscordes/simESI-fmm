#include "FileUtils.h"
#include "Parameters.h"

#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <sstream>

void deleteFile(const fs::path& fname) {
	if (!fs::exists(fname)) return;
	error_code ec;
	fs::remove(fname, ec);
	if (ec)throw runtime_error("Failed to delete file " + fname.string() + ": " + ec.message());
}

void createDir(const fs::path& dest) {
    error_code ec;
	fs::create_directory(dest, ec);
    if (ec)throw runtime_error("Failed to create directory " + dest.string() + ": " + ec.message());
}

void copyFile(const fs::path& source, const fs::path& dest) {
	if (!fs::exists(source)) throw runtime_error("No file " + source.string() + " to copy.");
	error_code ec;
	fs::copy_file(source, dest, fs::copy_options::overwrite_existing, ec);
	if (ec) throw runtime_error("Error copying file from " + source.string() + " to " + dest.string() + ": " + ec.message());
}

void copyFiles(const fs::path& idir, const fs::path& odir) {
	if (!fs::exists(idir)) throw runtime_error("No dir " + idir.string() + " to copy files from.");
	if (!fs::exists(odir)) throw runtime_error("No dir " + odir.string() + " to copy files to.");
	for (const auto& entry : fs::directory_iterator(idir)) {
		const fs::path& source = entry.path();
		const fs::path  dest = odir / source.filename();
		fs::copy_file(source, dest, fs::copy_options::skip_existing);
	}
}

void copyRecursive(const fs::path& source, const fs::path& dest) {
	if (!fs::exists(source)) throw runtime_error("No dir " + source.string() + " to copy.");
	error_code ec;
	fs::copy(source, dest, fs::copy_options::overwrite_existing | fs::copy_options::recursive, ec);
	if (ec) throw runtime_error("Error copying dir " + source.string() + " to " + dest.string() + ": " + ec.message());
}

void populateOutputDirs() {
    createDir("outputFiles"); // Create outputFiles if doesnt exist

    Parameters& params = Parameters::Get();
	string name;
	if (params.getSettingsPresent().test(OUTPUTDIRNAME))
		name = params.getOutputDirName();
	else {
		name = fs::path(params.getPDB()).stem().string();
		params.getOutputDirName() = name;
		params.getSettingsPresent().set(OUTPUTDIRNAME);
	}

    fs::path outputDir;
	if (!params.isCont()) {
		int count = 1;
		while (true) {
			outputDir = fs::path("outputFiles") / (name + '_' + to_string(count));
			if (!fs::exists(outputDir)) {
				createDir(outputDir);
				params.getOutputDir() = outputDir.string();

				params.getDataDir() = (outputDir / "data").string();
				createDir(params.getDataDir());

				params.getSimulationDir() = (outputDir / "simulation").string();
				createDir(params.getSimulationDir());

				break;
			}
			count++;
		}
	}
	else {
		outputDir = params.getDirCont();
		params.getOutputDir() = outputDir;
		params.getDataDir() = (outputDir / "data").string();
		params.getSimulationDir() = (outputDir / "simulation").string();
	}

	// Begin populating simulation dir
	fs::path simulationDir = params.getSimulationDir();
	fs::path resourceDir = "resources";
	fs::path ffDir = resourceDir / "charmm36.ff";
	fs::path mdpDir = resourceDir / "mdp_files";
	fs::path topDir = resourceDir / "top_files";

	copyFile(fs::path("coordinateFiles") / params.getPDB(), simulationDir / params.getPDB()); // .pdb

	if (params.getSettingsPresent().test(PKAPDB))
		copyFile(fs::path("coordinateFiles") / params.getPDB(), simulationDir / "pka.pdb"); // pka.pdb

	createDir(simulationDir / "charmm36.ff"); // FF files
	copyRecursive(ffDir, simulationDir / "charmm36.ff");

	// .mdp files, consider whether FMM used or not
	copyFile(mdpDir / "em.mdp", simulationDir / "em.mdp");
	if (params.isFMM()) {
		copyFile(mdpDir / "nvt_FMM.mdp", simulationDir / "nvt.mdp");
		copyFile(mdpDir / "prodrun_FMM.mdp", simulationDir / "prodrun.mdp");
	}
	else {
		copyFile(mdpDir / "nvt.mdp", simulationDir / "nvt.mdp");
		copyFile(mdpDir / "prodrun.mdp", simulationDir / "prodrun.mdp");
	}

	copyFiles(topDir, simulationDir); // .top filess
}

// Renames files to account for GROMACS checkpoint naming
void renameCheckpointFile(const string& ftype, int step, int num_restarts) {
	ostringstream oss;
	oss << setw(4) << setfill('0') << num_restarts;
	string s = to_string(step) + ".part" + oss.str() + ftype;
	string d = to_string(step) + ftype;

	try {
		if (fs::exists(s))
			fs::rename(s, d);
	}
	catch (const fs::filesystem_error& e) {
		throw runtime_error("renameCheckpointFile: Error when renaming file." + string(e.what()));
	}
}