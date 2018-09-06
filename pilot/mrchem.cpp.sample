/** The MRChem sandbox */

#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "getkw/Getkw.hpp"

#include "mrenv.h"
#include "mrchem.h"
#include "parallel.h"

Getkw mrchem::Input;
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace mrcpp;
using namespace mrchem;

int main(int argc, char **argv) {
    mpi::initialize(argc, argv);
    mrenv::initialize(argc, argv);

    Timer timer;

    // Do your stuff here
    
    timer.stop();
    double wt = timer.getWallTime();

    mrenv::finalize(wt);
    mpi::finalize();
    return 0;
}

