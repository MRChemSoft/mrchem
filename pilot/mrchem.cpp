/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

/** The MRChem sandbox */

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

int main(int argc, char **argv) {
    mpi::initialize();
    const auto json_inp = mrenv::fetch_json(argc, argv);
    mrenv::initialize(json_inp);

    Timer timer;

    // Do your stuff here
    println(0, json_inp.dump(2));

    timer.stop();
    double wt = timer.elapsed();

    mrenv::finalize(wt);
    mpi::finalize();
    return EXIT_SUCCESS;
}
