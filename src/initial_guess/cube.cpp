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

#include "cube.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "parallel.h"
#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <highfive/H5Easy.hpp>
#include <nlohmann/json.hpp>

#include "analyticfunctions/CUBEfunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using nlohmann::json;

namespace mrchem {

namespace initial_guess {
namespace cube {
bool setup(OrbitalVector &Phi, double prec, const std::string &file_p, const std::string &file_a, const std::string &file_b) {
    if (Phi.size() == 0) return false;

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation   ", "Compute initial orbitals");
    print_utils::text(0, "Method        ", "Project cube file molecular orbitals");
    print_utils::text(0, "Precision     ", print_utils::dbl_to_str(prec, 5, true));
    if (orbital::size_singly(Phi)) {
        print_utils::text(0, "Restricted    ", "False");
        print_utils::text(0, "MO alpha file ", file_a);
        print_utils::text(0, "MO beta file  ", file_b);
    } else {
        print_utils::text(0, "Restricted    ", "True");
        print_utils::text(0, "MO file ", file_p);
    }
    mrcpp::print::separator(0, '~', 2);

    // Separate alpha/beta from paired orbitals
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);

    // Project paired, alpha and beta separately
    auto success = true;
    success &= initial_guess::cube::project_mo(Phi, prec, file_p);
    success &= initial_guess::cube::project_mo(Phi_a, prec, file_a);
    success &= initial_guess::cube::project_mo(Phi_b, prec, file_b);

    // Collect orbitals into one vector
    for (auto &phi_a : Phi_a) Phi.push_back(phi_a);
    for (auto &phi_b : Phi_b) Phi.push_back(phi_b);

    return success;
}

bool project_mo(OrbitalVector &Phi, double prec, const std::string &mo_file) {
    if (Phi.size() == 0) return true;

    Timer t_tot;
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 2;
    auto w1 = 5;
    auto w2 = w0 * 2 / 9;
    auto w3 = w0 - w1 - 3 * w2;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w3) << "Norm";
    o_head << std::setw(w2 + 1) << "Nodes";
    o_head << std::setw(w2) << "Size";
    o_head << std::setw(w2) << "Time";

    mrcpp::print::header(1, "CUBE Initial Guess");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    json cube_inp;
    std::ifstream ifs(mo_file, std::ios_base::in);
    ifs >> cube_inp;
    ifs.close();

    auto CUBEVector = getCUBEFunction(cube_inp);

    bool success = true;
    for (int i = 0; i < Phi.size(); i++) {
        Timer t_i;
        if (mpi::my_orb(Phi[i])) {
            CUBEfunction phi_i = CUBEVector[i];
            Phi[i].alloc(NUMBER::Real);
            mrcpp::project(prec, Phi[i].real(), phi_i);
            std::stringstream o_txt;
            o_txt << std::setw(w1 - 1) << i;
            o_txt << std::setw(w3) << print_utils::dbl_to_str(Phi[i].norm(), pprec, true);
            print_utils::qmfunction(1, o_txt.str(), Phi[i], t_i);
        }
    }
    mpi::barrier(mpi::comm_orb);
    mrcpp::print::footer(1, t_tot, 2);
    return success;
}

std::vector<mrchem::CUBEfunction> getCUBEFunction(const json &json_cube) {
    std::vector<mrchem::CUBEfunction> CUBEVector;
    for (const auto &item : json_cube.items()) {
        auto Header = item.value()["Header"];

        auto origin = Header["origin"];
        auto N_steps = Header["N_steps"];
        auto Voxel_axes = Header["Voxel_axes"];

        for (const auto &x : item.value()["CUBE_data"].items()) {
            // the data is saved as a vector of vectors indexing as
            // CUBE_data[ID][x_val*n_steps[1]*n_steps[2] + y_val*n_steps[2] + z_val]
            CUBEVector.emplace_back(N_steps, origin, Voxel_axes, x.value());
        }
    }

    return CUBEVector;
}
} // namespace cube

namespace hdf5 {
bool setup(OrbitalVector &Phi, double prec, const std::string &h5fname) {
    if (Phi.size() == 0) return false;

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation   ", "Compute initial orbitals");
    print_utils::text(0, "Method        ", "Project HDF5 molecular orbitals");
    print_utils::text(0, "Precision     ", print_utils::dbl_to_str(prec, 5, true));
    auto is_restricted = orbital::size_singly(Phi) ? "False" : "True";
    // we have one HDF5 file regardless of restricted or unrestricted
    print_utils::text(0, "Restricted    ", is_restricted);
    print_utils::text(0, "MO file ", h5fname);
    mrcpp::print::separator(0, '~', 2);

    // load HDF5 file
    H5Easy::File h5file(h5fname, H5Easy::File::ReadOnly);

    // Separate alpha/beta from paired orbitals
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);

    // Project paired, alpha and beta separately
    auto success = true;
    success &= hdf5::project_mo(Phi, prec, h5file, "p");
    success &= hdf5::project_mo(Phi_a, prec, h5file, "a");
    success &= hdf5::project_mo(Phi_b, prec, h5file, "b");

    // Collect orbitals into one vector
    for (auto &phi_a : Phi_a) Phi.push_back(phi_a);
    for (auto &phi_b : Phi_b) Phi.push_back(phi_b);

    return success;
}

bool project_mo(OrbitalVector &Phi, double prec, const H5Easy::File & h5file, const std::string& spin_lbl) {
    if (Phi.size() == 0) return true;

    Timer t_tot;
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 2;
    auto w1 = 5;
    auto w2 = w0 * 2 / 9;
    auto w3 = w0 - w1 - 3 * w2;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w3) << "Norm";
    o_head << std::setw(w2 + 1) << "Nodes";
    o_head << std::setw(w2) << "Size";
    o_head << std::setw(w2) << "Time";

    mrcpp::print::header(1, "HDF5 Initial Guess");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    auto CUBEVector = getCUBEFunction(h5file, spin_lbl);

    bool success = true;
    for (int i = 0; i < Phi.size(); i++) {
        Timer t_i;
        if (mpi::my_orb(Phi[i])) {
            CUBEfunction phi_i = CUBEVector[i];
            Phi[i].alloc(NUMBER::Real);
            mrcpp::project(prec, Phi[i].real(), phi_i);
            std::stringstream o_txt;
            o_txt << std::setw(w1 - 1) << i;
            o_txt << std::setw(w3) << print_utils::dbl_to_str(Phi[i].norm(), pprec, true);
            print_utils::qmfunction(1, o_txt.str(), Phi[i], t_i);
        }
    }
    mpi::barrier(mpi::comm_orb);
    mrcpp::print::footer(1, t_tot, 2);
    return success;
}

std::vector<mrchem::CUBEfunction> getCUBEFunction(const H5Easy::File &h5file, const std::string & spin_lbl) {
    std::vector<mrchem::CUBEfunction> CUBEVector;

    // attributes
    Eigen::Vector3i num_points;
    Eigen::Vector3d origin;
    Eigen::Matrix3d step_size;
    // data
    Eigen::VectorXd orbital;

    for (auto i = 0; i < h5file.getGroup("ground_state").getNumberObjects(); ++i) {

        std::ostringstream lbl;
        lbl << "/ground_state/" << std::setfill('0') << std::setw(4) << i << "/mo_" << spin_lbl;

        // load attributes for each MO
        num_points = H5Easy::loadAttribute<Eigen::Vector3i>(h5file, lbl.str(), "num_points");
        origin = H5Easy::loadAttribute<Eigen::Vector3d>(h5file, lbl.str(), "origin");
        step_size = H5Easy::loadAttribute<Eigen::Matrix3d>(h5file, lbl.str(), "step_size");

        orbital = H5Easy::load<Eigen::VectorXd>(h5file, lbl.str());

        CUBEVector.emplace_back(num_points, origin, step_size, orbital);
    }

    return CUBEVector;
}
} // namespace hdf5
} // namespace initial_guess
} // namespace mrchem
