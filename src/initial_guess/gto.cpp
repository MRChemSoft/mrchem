#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "Getkw.h"

#include "mrenv.h"
#include "mrchem.h"
#include "parallel.h"

#include "Molecule.h"
#include "Orbital.h"
#include "gto_guess.h"

Getkw mrchem::Input;
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace mrcpp;
using namespace mrchem;

/** @file gto.cpp
 *
 * Standalone executable (gto-guess) for reading a GTO initial guess and
 * writing the resulting MW orbitals to disk.
 *
 * Requires the following input files (file names can be changed in input):
 * @mrchem.inp: regular input file, parsed through getkw (./mrchem -D)
 * initial_guess/mrchem.bas: basis set file (LSDalton format)
 * initial_guess/mrchem.moa: MO matrix (alpha orbitals or closed-shell)
 * initial_guess/mrchem.mob: MO matrix (beta orbitals)
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
    std::string bas_file = Input.get<string>("Files.basis_set");
    std::string moa_file = Input.get<string>("Files.mo_mat_a");
    std::string mob_file = Input.get<string>("Files.mo_mat_b");

    if (not (scf_guess == "GTO" or scf_guess == "gto" or scf_guess == "none"))
        MSG_FATAL("Invalid initial guess");

    // Setting up molecule
    Molecule mol(mol_coords, mol_charge, mol_multiplicity);
    mol.printGeometry();

    // Setting up orbitals
    OrbitalVector Phi;
    if (wf_restricted) {
        Phi = gto_guess::initial_guess(prec, mol, bas_file, moa_file);
    } else {
        Phi = gto_guess::initial_guess(prec, mol, bas_file, moa_file, mob_file);
    }
    orbital::save_orbitals(Phi, orb_file);
    orbital::free(Phi);

    timer.stop();
    mrenv::finalize(timer.getWallTime());
    mpi::finalize();

    return 0;
}

