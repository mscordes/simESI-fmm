#pragma once
#include "Writer.h"
#include "State.h"

class ChargeWriter : public Writer {
    public:
        void writeHeader();
        void write(const State& state);
};