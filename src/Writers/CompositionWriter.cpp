#include "CompositionWriter.h"
#include "Constants.h"
#include "Parameters.h"

void CompositionWriter::writeHeader() 
{
    outFile << "# ";
    for (const string& resType : Constants::topOrder)
        outFile << setw(6) << resType; 
    outFile << setw(6) << "gSOL";
    outFile << endl;
}

void CompositionWriter::write(
    const State& state, 
    const unordered_set<shared_ptr<Residue>>& gWaters) 
{
    static const bool humid = Parameters::Get().isHumid();
    outFile << "  ";
    for (const string& resType : Constants::topOrder) {
        auto it = state.residueSet.find(resType);
        if (it == state.residueSet.end()) outFile << setw(6) << "0";
        else {
            if (humid && resType == "SOL") 
                outFile << setw(6) << it->second.size() - gWaters.size();
            else 
                outFile << setw(6) << it->second.size();
        }
    }

    if (humid)
        outFile << setw(6) << gWaters.size();
    else
        outFile << setw(6) << "0";

    outFile << endl;
}