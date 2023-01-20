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

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "mrchem.h"
#include "parallel.h"

#include "utils/print_utils.h"

#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"

mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

void initialize_mra();
void finalize_mra();

int main(int argc, char *argv[]) {
    // global setup
    initialize_mra();
    mrchem::mpi::bank_size = 1;            // Will be set to zero if world_size == 1
    mrchem::mpi::numerically_exact = true; // Required for MPI invariant results
    mrchem::mpi::initialize();

    // print MPI info
    mrcpp::Printer::init(0, mrchem::mpi::world_rank, mrchem::mpi::world_size);
    mrcpp::print::separator(0, '-');
    mrchem::print_utils::scalar(0, "MPI processes  ", mrchem::mpi::world_size, "(" + std::to_string(mrchem::mpi::bank_size) + " bank)", 0, false);
    mrchem::print_utils::scalar(0, "OpenMP threads ", mrchem::omp::n_threads, "", 0, false);
    mrchem::print_utils::scalar(0, "Total cores    ", mrchem::mpi::world_size * mrchem::omp::n_threads, "", 0, false);
    mrcpp::print::separator(0, '-');

    // run tests
    int result = Catch::Session().run(argc, argv);

    // global cleanup
    finalize_mra();
    mrchem::mpi::finalize();
    return (result < 0xff ? result : 0xff);
}

void initialize_mra() {
    // Defining resolution levels 2^{-n}
    int min_scale = -5;
    int max_depth = 25;

    // Constructing world box
    std::array<int, 3> boxes{2, 2, 2};
    std::array<int, 3> corner{-1, -1, -1};
    std::array<double, 3> sfac{1.0, 1.0, 1.0};
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes, sfac);

    // Constructing basis
    int order = 5;
    mrcpp::InterpolatingBasis basis(order);

    // Constructing global MRA
    mrchem::MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
}

void finalize_mra() {
    if (mrchem::MRA != nullptr) delete mrchem::MRA;
    mrchem::MRA = nullptr;
}
