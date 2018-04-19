#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "Getkw.h"

#include "mrenv.h"
#include "mrchem.h"
#include "parallel.h"
#include "hydrogen_guess.h"

#include "HydrogenFunction.h"
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

    int ig_zeta = 0;
         if (scf_guess == "SZ") { ig_zeta = 1; }
    else if (scf_guess == "DZ") { ig_zeta = 2; }
    else if (scf_guess == "TZ") { ig_zeta = 3; }
    else if (scf_guess == "QZ") { ig_zeta = 4; }
    else { MSG_FATAL("Invalid initial guess"); }

    // Setting up molecule
    Molecule mol(mol_coords, mol_charge, mol_multiplicity);
    mol.printGeometry();

    // Setting up orbitals
    OrbitalVector Phi = hydrogen_guess::initial_guess(prec, mol, wf_restricted, ig_zeta);
    orbital::save_orbitals(Phi, orb_file);
    orbital::free(Phi);

    timer.stop();
    mrenv::finalize(timer.getWallTime());
    mpi::finalize();

    return 0;
}

