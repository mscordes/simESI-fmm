#include "TemperatureWriter.h"

void TemperatureWriter::writeHeader() {
    outFile << "# Temperature (K)" << endl;
}

void TemperatureWriter::write(double temperature) {
    outFile << temperature << endl;
}