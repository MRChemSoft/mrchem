/** The MRChem sandbox */

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"

#include "analyticfunctions/HydrogenFunction.h"
#include "chemistry/Cavity.h"
#include "chemistry/Element.h"
#include "chemistry/Nucleus.h"
#include "chemistry/PeriodicTable.h"
#include "chemistry/Permittivity.h"
#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/two_electron/ReactionOperator.h"
#include "qmoperators/two_electron/SCRF.h"

// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

int main(int argc, char **argv) {
    const auto json_inp = mrenv::fetch_json(argc, argv);
    mrenv::initialize(json_inp);

    Timer timer;

    // Do your stuff here
    println(0, json_inp.dump(2));

    timer.stop();
    double wt = timer.elapsed();

    mrenv::finalize(wt);
    mpi::finalize();
    return EXIT_SUCCESS;
}
