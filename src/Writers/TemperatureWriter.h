#pragma once
#include "Writer.h"
#include "State.h"

class TemperatureWriter : public Writer {
    public:
        void writeHeader();
        void write(double temperature);
};