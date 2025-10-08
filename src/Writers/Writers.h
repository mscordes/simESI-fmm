#pragma once
#include "ChargeWriter.h"
#include "CompositionWriter.h"
#include "ExchangeWriter.h"
#include "TemperatureWriter.h"
#include "TimerWriter.h"
#include "State.h"

struct Writers {
    fs::path            dataDir;
    ChargeWriter        chargeWriter;
    CompositionWriter   compositionWriter;
    ExchangeWriter      exchangeWriter;
    TemperatureWriter   temperatureWriter;
    TimerWriter         timerWriter;

    Writers(const fs::path& dataDir, const State& state);
    ~Writers() = default;
};