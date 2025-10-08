#include "Writers.h"  
#include "Parameters.h"

Writers::Writers(const fs::path& dataDir_, const State& state)
    : dataDir(dataDir_),
    chargeWriter(dataDir_ / "charge.txt"),
    compositionWriter(dataDir_ / "composition.txt"),
    exchangeWriter(dataDir_ / "exchanges.txt"),
    temperatureWriter(dataDir_ / "temperature.txt"),
    timerWriter(dataDir_ / "timer.txt")
{
    if (!Parameters::Get().isCont()) {
        chargeWriter.writeHeader();
        compositionWriter.writeHeader();
        temperatureWriter.writeHeader();
        timerWriter.writeHeader();

        chargeWriter.write(state);
        compositionWriter.write(state, getGaseousWaters(state));
    }
}