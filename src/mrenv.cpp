#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "getkw/Getkw.hpp"

#include "mrenv.h"
#include "mrchem.h"
#include "parallel.h"

using namespace std;
using namespace mrcpp;
using namespace mrchem;

namespace mrenv {

void initialize(int argc, char **argv) {
    const char *infile = 0;
    if (argc == 1) {
        infile = "STDIN";
    } else if (argc == 2) {
        infile = argv[1];
    } else {
        MSG_ERROR("Ivalid number of arguments!");
    }

    // Initialize input
    Input = Getkw(infile, false, true);

    // Initialize printing
    int printlevel = Input.get<int>("printlevel");
    bool teletype = Input.get<bool>("teletype");
    if (teletype) {
        Printer::init(printlevel, mpi::orb_rank, mpi::orb_size, "mrchem");
    } else {
        Printer::init(printlevel, mpi::orb_rank, mpi::orb_size);
    }
    Printer::setPrecision(15);

    // Initialize world box
    int min_scale = Input.get<int>("MRA.min_scale");
    int max_scale = Input.get<int>("MRA.max_scale");
    vector<int> corner = Input.getIntVec("MRA.corner");
    vector<int> boxes = Input.getIntVec("MRA.boxes");
    BoundingBox<3> world(min_scale, corner.data(), boxes.data());

    // Initialize scaling basis
    int order = Input.get<int>("MRA.order");
    string btype = Input.get<string>("MRA.basis_type");

    int max_depth = max_scale - min_scale;
    if (min_scale < MinScale) MSG_FATAL("Root scale too large");
    if (max_scale > MaxScale) MSG_FATAL("Max scale too large");
    if (max_depth > MaxDepth) MSG_FATAL("Max depth too large");

    // Initialize global MRA
    if (btype == "I") {
        InterpolatingBasis basis(order);
        MRA = new MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else if (btype == "L") {
        LegendreBasis basis(order);
        MRA = new MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else {
        MSG_FATAL("Invalid basis type!");
    }

    println(0,endl << endl);
    println(0,"************************************************************");
    println(0,"***     __  __ ____   ____ _                             ***");
    println(0,"***    |  \\/  |  _ \\ / ___| |__   ___ _ __ ___           ***");
    println(0,"***    | |\\/| | |_) | |   | '_ \\ / _ \\ '_ ` _ \\          ***");
    println(0,"***    | |  | |  _ <| |___| | | |  __/ | | | | |         ***");
    println(0,"***    |_|  |_|_| \\_\\\\____|_| |_|\\___|_| |_| |_|         ***");
    println(0,"***                                                      ***");
    println(0,"***    VERSION " << PROGRAM_VERSION << " (rev. " << GIT_REVISION << ")                      ***");
    println(0,"***                                                      ***");
    println(0,"***    Stig Rune Jensen <stig.r.jensen@uit.no>           ***");
    println(0,"***    Luca Frediani    <luca.frediani@uit.no>           ***");
    println(0,"***    Peter Wind       <peter.wind@uit.no>              ***");
    println(0,"***                                                      ***");
    println(0,"************************************************************");
    println(0,endl);

    if (mpi::orb_size > 1 or omp::n_threads > 1) {
        println(0,"+++ Parallel execution: ");
        println(0,"  MPI hosts available     : " << mpi::orb_size);
        println(0,"  Threads/host            : " << omp::n_threads);
        println(0,"  Total used CPUs         : " << mpi::orb_size*omp::n_threads);
        println(0,"");
    } else {
        println(0,"+++ Serial execution" << endl);
    }

    Printer::printEnvironment();

    int oldprec = Printer::setPrecision(6);
    MRA->print();
    Printer::setPrecision(oldprec);
}

void finalize(double wt) {
    // Delete global MRA
    if (MRA != 0) delete MRA;
    MRA = 0;

    Printer::setPrecision(6);
    println(0,endl);
    println(0,"************************************************************");
    println(0,"***                                                      ***");
    println(0,"***                    Exiting MRChem                    ***");
    println(0,"***                                                      ***");
    println(0,"***               Wall time: " << wt << "                ***");
    println(0,"***                                                      ***");
    println(0,"************************************************************");
    println(0,endl);
}

} //namespace mrenv
