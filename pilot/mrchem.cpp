/** The MRChem sandbox */

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

int main(int argc, char **argv) {
    mpi::initialize();
    const auto json_input = mrenv::fetch_json(argc, argv);
    mrenv::initialize(json_input);

    Timer timer;

    // Do your stuff here
    println(0, json_input.dump(2));

    timer.stop();
    double wt = timer.elapsed();

    mrenv::finalize(wt);
    mpi::finalize();
    return 0;
}
