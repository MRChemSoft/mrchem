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
    const auto json_inp = mrenv::fetch_json(argc, argv);

    json json_out;
    mrenv::initialize(json_inp);
    const auto &mol_inp = json_inp["molecule"];
    const auto &scf_inp = json_inp["scf_calculation"];
    const auto &rsps_inp = json_inp["rsp_calculations"];

    Timer timer;
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    json_out["scf_calculation"] = driver::scf::run(scf_inp, mol);
    if (json_out["scf_calculation"]["success"]) {
        for (const auto &rsp_inp : rsps_inp) driver::rsp::run(rsp_inp, mol);
    }
    mpi::barrier(mpi::comm_orb);
    driver::print_properties(mol);
    mrenv::finalize(timer.elapsed());
    mrenv::dump_json(json_inp, json_out);

    mpi::finalize();
    return EXIT_SUCCESS;
}
