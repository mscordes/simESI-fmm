#pragma once
#include "Writer.h"
#include "State.h"

class CompositionWriter : public Writer {
    public:
        void writeHeader();
        void write(
            const State& state,
            const unordered_set<shared_ptr<Residue>>& gWaters);
};