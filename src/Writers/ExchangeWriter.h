#pragma once
#include "Writer.h"
#include "Exchange.h"
#include "State.h"

#include <vector>

class ExchangeWriter : public Writer {
    public:
        void write(const unordered_set<shared_ptr<Exchange>>& exchanges);
};