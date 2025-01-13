#include "Core.h"
#include <filesystem>

int main(int argc, char** argv) {

    // Parse user arguments
	Core::Config config = Core::getArgs(argc, argv);

    // Make and populate trial dirs
    std::filesystem::path trialPath = Core::makeDirs(config);

    // Run simulation
	Core::simulation(config, trialPath);

    return 0;
}