/** \mainpage The MRChem main program
 *
 * \author Stig Rune Jensen
 *
 * \version 1.0
 *
 * \par Copyright:
 * GPLv4
 *
 */

#include "Getkw.h"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

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

    timer.stop();
    double wt = timer.getWallTime();

    mrenv::finalize(wt);
    mpi::finalize();
    return 0;
}

