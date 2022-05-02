/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
#include <fstream>
#include <string>

#include "analyticfunctions/CUBEfunction.h"
#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

std::shared_ptr<mrchem::CUBEfunction> getCUBEFunction(const json &json_cube, int i);
int main(int argc, char **argv) {
    mpi::initialize();
    const auto json_inp = mrenv::fetch_json(argc, argv);
    mrenv::initialize(json_inp);

    Timer timer;

    // Do your stuff here

    auto prec = 1.0e-3;

    std::string mo_file = "pilot_vec_p.cube";
    json cube_inp;
    std::ifstream ifs(mo_file, std::ios_base::in);
    ifs >> cube_inp;
    ifs.close();
    Orbital phi(1.0);
    phi.alloc(NUMBER::Real);
    for (int i = 0; i < 10; i++) {
        auto phi_i = getCUBEFunction(cube_inp, i % 5);
        mrcpp::project(prec, phi.real(), *phi_i);
        std::cout << "first norm of phi at cycle " << i << ": " << phi.norm() << "\n";
        phi.free(NUMBER::Real);

        phi.alloc(NUMBER::Real);
        mrcpp::project(prec, phi.real(), *phi_i);
        std::cout << "first norm of phi at cycle " << i << ": " << phi.norm() << "\n";
        phi.free(NUMBER::Real);

        phi.alloc(NUMBER::Real);
        mrcpp::project(prec, phi.real(), *phi_i);
        std::cout << "first norm of phi at cycle " << i << ": " << phi.norm() << "\n";
        phi.free(NUMBER::Real);

        phi.alloc(NUMBER::Real);
        mrcpp::project(prec, phi.real(), *phi_i);
        std::cout << "first norm of phi at cycle " << i << ": " << phi.norm() << "\n";

        OrbitalVector orb_vec;

        for (auto j = 0; j < 5; j++) {
            Orbital psi(1.0);
            psi.alloc(NUMBER::Real);
            auto psi_j = getCUBEFunction(cube_inp, j % 5);
            mrcpp::project(prec, psi.real(), *psi_j);
            std::cout << "first norm of psi at cycle " << i << " and " << j << ": " << psi.norm() << "\n";
            orb_vec.push_back(psi);
        }
        std::cout << "norms of the orbitalvector: " << orbital::get_norms(orb_vec) << "\n";
    }
    // Did my stuff there

    println(0, json_inp.dump(2));

    timer.stop();
    double wt = timer.elapsed();

    mrenv::finalize(wt);
    mpi::finalize();
    return EXIT_SUCCESS;
}

std::shared_ptr<mrchem::CUBEfunction> getCUBEFunction(const json &json_cube, int index) {
    auto Header = json_cube[index]["Header"];
    auto N_atoms = Header["N_atoms"];
    auto origin = Header["origin"];
    auto N_steps = Header["N_steps"];
    auto Voxel_axes = Header["Voxel_axes"];
    auto Z_n = Header["Z_n"];
    auto atom_charges = Header["atom_charges"];
    auto atom_coords = Header["atom_coords"];
    auto N_vals = Header["N_vals"];
    auto CUBE_data = json_cube[index]["CUBE_data"][0];

    return std::make_shared<CUBEfunction>(N_atoms, N_vals, N_steps, origin, Voxel_axes, Z_n, CUBE_data, atom_charges, atom_coords);
}
