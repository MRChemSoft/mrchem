/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "getkw/Getkw.hpp"

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"
#include "sad.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

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

// clang-format off
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
    std::string scf_guess = Input.get<std::string>("SCF.initial_guess");
    std::string orb_file = Input.get<std::string>("Files.start_orbitals");

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

    timer.stop();
    mrenv::finalize(timer.getWallTime());
    mpi::finalize();

    return 0;
}
// clang-format on
