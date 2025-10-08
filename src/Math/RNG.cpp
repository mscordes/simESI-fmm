#include "RNG.h"

namespace RNG {
    mt19937 gen{ random_device{}() };
}