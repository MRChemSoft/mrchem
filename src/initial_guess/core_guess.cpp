/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "chemistry/Molecule.h"
#include "core.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

Getkw mrchem::input;
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace mrcpp;
using namespace mrchem;

/** @file core_guess.cpp
 *
 * Standalone executable (core-guess) for setting up an initial guess
 * from hydrogen orbitals and writing the resulting MW orbitals to disk.
 *
 * Sets up an AO basis of hydrogen functions with the given zeta quality
 * (SZ, DZ, TZ, QZ), computes and diagonalizes the core Hamiltonian matrix,
 * and fills the resulting orbitals by the Aufbau principle.
 *
 * Requires the following input files:
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
    double prec = input.get<double>("rel_prec");
    bool wf_restricted = input.get<bool>("wavefunction.restricted");
    int mol_charge = input.get<int>("molecule.charge");
    int mol_multiplicity = input.get<int>("molecule.multiplicity");
    std::vector<std::string> mol_coords = input.getData("molecule.coords");
    std::string scf_guess = input.get<std::string>("scf.initial_guess");
    std::string orb_file = input.get<std::string>("files.start_orbitals");

    int ig_zeta = 0;
         if (scf_guess == "core_sz") { ig_zeta = 1; }
    else if (scf_guess == "core_dz") { ig_zeta = 2; }
    else if (scf_guess == "core_tz") { ig_zeta = 3; }
    else if (scf_guess == "core_qz") { ig_zeta = 4; }
    else { MSG_FATAL("Invalid initial guess"); }

    // Setting up molecule
    Molecule mol(mol_coords, mol_charge, mol_multiplicity);
    mol.printGeometry();

    // Setting up orbitals
    OrbitalVector Phi = initial_guess::core::setup(prec, mol, wf_restricted, ig_zeta);
    orbital::save_orbitals(Phi, orb_file);

    timer.stop();
    mrenv::finalize(timer.getWallTime());
    mpi::finalize();

    return 0;
}
// clang-format on