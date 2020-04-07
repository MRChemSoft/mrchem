/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <MRCPP/Timer>

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
    if (driver::scf::run(json_scf, mol)) {
        for (const auto &json_rsp : json_rsps) driver::rsp::run(json_rsp, mol);
    }
    driver::print_properties(mol);
    mpi::barrier(mpi::comm_orb);
    mrenv::finalize(timer.elapsed());

    mpi::finalize();
    return EXIT_SUCCESS;
}
