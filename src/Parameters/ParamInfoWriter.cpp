#include "ParamInfoWriter.h"
#include "Parameters.h"
#include <stdexcept>

using namespace std;

ParamInfoWriter::ParamInfoWriter(const fs::path& fileName_)
    : fileName(fileName_) {
}

void ParamInfoWriter::write() const {
    ofstream outFile(fileName);
    if (!outFile) {
        throw runtime_error("Could not open output file: " + fileName.string());
    }
    Parameters::Get().print(outFile);
}