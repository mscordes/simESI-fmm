#pragma once
#include "Writer.h"
#include "State.h"
#include "Timer.h"

class TimerWriter : public Writer {
    public:
        void writeHeader();
        void write(
            const chrono::steady_clock::time_point& tInit,
            const chrono::steady_clock::time_point& tStart,
            const chrono::steady_clock::time_point& tExchange,
            const chrono::steady_clock::time_point& tTop,
            const chrono::steady_clock::time_point& tGrompp,
            const chrono::steady_clock::time_point& tMD);
};