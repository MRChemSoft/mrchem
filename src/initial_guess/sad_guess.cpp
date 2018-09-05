#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "Getkw.h"

#include "mrenv.h"
#include "mrchem.h"
#include "parallel.h"
#include "initial_guess/sad.h"

#include "Molecule.h"
#include "Orbital.h"
#include "orbital_utils.h"

Getkw mrchem::Input;
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace mrcpp;
using namespace mrchem;

/** @file sad_guess.cpp
 *
 * Standalone executable (sad-guess) for reading a GTO initial guess and
 * writing the resulting MW orbitals to disk.
 *
 * Requires the following input files (file names can be changed in input):
 * @mrchem.inp: regular input file, parsed through getkw (./mrchem -D)
 *
 * Produces the following output files (file names can be changed in input):
 * orbitals/phi_0.meta: orbital meta data
 * orbitals/phi_0_re.tree: MW representation of real part
 * orbitals/phi_1.meta: orbital meta data
 * orbitals/phi_1_re.tree: MW representation of real part
 */

int main(int argc, char **argv) {
    mpi::initialize(argc, argv);
    mrenv::initialize(argc, argv);
    Timer timer;

    // Reading input
    double prec = Input.get<double>("rel_prec");
    bool wf_restricted = Input.get<bool>("WaveFunction.restricted");
    int mol_charge = Input.get<int>("Molecule.charge");
    int mol_multiplicity = Input.get<int>("Molecule.multiplicity");
    std::vector<std::string> mol_coords = Input.getData("Molecule.coords");
    std::string scf_guess = Input.get<string>("SCF.initial_guess");
    std::string orb_file = Input.get<string>("Files.start_orbitals");

    int ig_zeta = 0;
         if (scf_guess == "SAD_SZ") { ig_zeta = 1; }
    else if (scf_guess == "SAD_DZ") { ig_zeta = 2; }
    else if (scf_guess == "SAD_TZ") { ig_zeta = 3; }
    else if (scf_guess == "SAD_QZ") { ig_zeta = 4; }
    else { MSG_FATAL("Invalid initial guess"); }

    // Setting up molecule
    Molecule mol(mol_coords, mol_charge, mol_multiplicity);
    mol.printGeometry();

    // Setting up orbitals
    OrbitalVector Phi = initial_guess::sad::setup(prec, mol, wf_restricted, ig_zeta);
    orbital::save_orbitals(Phi, orb_file);
    orbital::free(Phi);

    timer.stop();
    mrenv::finalize(timer.getWallTime());
    mpi::finalize();

    return 0;
}
