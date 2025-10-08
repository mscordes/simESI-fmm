#include "TimerWriter.h"

void TimerWriter::writeHeader() {
    outFile << "#   simESI(s)  pdb2gmx(s)   grompp(s)    mdrun(s)   total(hr)" << endl;
}

void TimerWriter::write(
    const chrono::steady_clock::time_point& tInit,
    const chrono::steady_clock::time_point& tStart,
    const chrono::steady_clock::time_point& tExchange,
    const chrono::steady_clock::time_point& tTop,
    const chrono::steady_clock::time_point& tGrompp,
    const chrono::steady_clock::time_point& tMD)
{
    outFile << fixed << setprecision(3) << right
        << setw(13) << seconds(tExchange - tStart)
        << setw(12) << seconds(tTop - tExchange)
        << setw(12) << seconds(tGrompp - tTop)
        << setw(12) << seconds(tMD - tGrompp)
        << setw(12) << duration_cast<chrono::duration<double>>(now() - tInit).count() / 3600.0
        << endl;
}