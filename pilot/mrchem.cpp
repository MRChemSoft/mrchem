/** The MRChem sandbox */

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

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
