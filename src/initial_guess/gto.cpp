#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "Getkw.h"

#include "mrenv.h"
#include "mrchem.h"
#include "parallel.h"
#include "gto_guess.h"

#include "Molecule.h"
#include "Orbital.h"

Getkw mrchem::Input;
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace mrcpp;
using namespace mrchem;

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

    if (not (scf_guess == "GTO" or scf_guess == "none")) MSG_FATAL("Invalid initial guess");

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

