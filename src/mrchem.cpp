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
#include "MRCPP/Timer"
#include "parallel.h"

#include "driver.h"
#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

#include "chemistry/Molecule.h"

// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

int main(int argc, char **argv) {
    const auto json_input = mrenv::fetch_input(argc, argv);

    mrenv::initialize(json_input);
    const auto &json_mol = json_input["molecule"];
    const auto &json_scf = json_input["scf_calculation"];
    const auto &json_rsps = json_input["rsp_calculations"];

    Timer timer;
    Molecule mol;
    driver::init_molecule(json_mol, mol);
    if (driver::run_scf(json_scf, mol)) {
        for (const auto &json_rsp : json_rsps) driver::run_rsp(json_rsp, mol);
    }
    driver::print_properties(mol);
    mpi::barrier(mpi::comm_orb);
    mrenv::finalize(timer.elapsed());

    mpi::finalize();
    return EXIT_SUCCESS;
}
