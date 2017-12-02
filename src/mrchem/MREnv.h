#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "parallel.h"
#include "constants.h"
#include "Timer.h"

class MREnv {
public:
    static void initializeMRCPP(int argc, char **argv);
    static void finalizeMRCPP(const Timer t);
    static void initializeMRA();
};

