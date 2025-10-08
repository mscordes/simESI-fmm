#include "Subprocess.h"
#include "Parameters.h"

#include <cstdio>  
#include <stdexcept>

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#define DEV_NULL "NUL" 
#else
#define POPEN popen
#define PCLOSE pclose
#define DEV_NULL "/dev/null"
#endif

// Call terminal with a given command and inputs
void subprocess(string& command, const vector<string>& inputs) {
    static bool verbose = Parameters::Get().isVerbose();

    if (!verbose) command += " > " DEV_NULL " 2>&1";

    FILE* pipe = POPEN(command.c_str(), "w");
    if (!pipe) throw runtime_error("Failed to run command: " + command);

    for (const auto& input : inputs) 
        fputs((input + "\n").c_str(), pipe);

    if (PCLOSE(pipe) == -1) throw runtime_error("Failed to close command pipe: " + command);
}