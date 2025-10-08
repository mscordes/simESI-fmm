#include "ChargeWriter.h"
#include "Charge.h"

void ChargeWriter::writeHeader() {
    outFile << "# Protein   System" << endl;
}

void ChargeWriter::write(const State& state) {
    double protCharge = getProteinCharge(state);
    double netCharge = getNetCharge(state, protCharge);
    outFile << fixed << setprecision(3)
        << setw(9) << protCharge
        << setw(9) << netCharge
        << endl;
}