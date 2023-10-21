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

#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "mw.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace mw {

bool project_mo(OrbitalVector &Phi, double prec, const std::string &mo_file);

} // namespace mw
} // namespace initial_guess

bool initial_guess::mw::setup(OrbitalVector &Phi, double prec, const std::string &file_p, const std::string &file_a, const std::string &file_b) {
    if (Phi.size() == 0) return false;

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation   ", "Compute initial orbitals");
    print_utils::text(0, "Method        ", "Project MW molecular orbitals");
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
    success &= initial_guess::mw::project_mo(Phi, prec, file_p);
    success &= initial_guess::mw::project_mo(Phi_a, prec, file_a);
    success &= initial_guess::mw::project_mo(Phi_b, prec, file_b);

    // Collect orbitals into one vector
    Phi = orbital::adjoin(Phi, Phi_a);
    Phi = orbital::adjoin(Phi, Phi_b);

    return success;
}

bool initial_guess::mw::project_mo(OrbitalVector &Phi, double prec, const std::string &mo_file) {
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

    mrcpp::print::header(1, "MW Initial Guess");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    bool success = true;
    for (int i = 0; i < Phi.size(); i++) {
        Timer t_i;
        if (mrcpp::mpi::my_orb(Phi[i])) {
            std::stringstream orbname;
            orbname << mo_file << "_idx_" << i;

            Orbital phi_i;
            orbital::loadOrbital(orbname.str(), phi_i);
            if (phi_i.squaredNorm() < 0.0) {
                MSG_ERROR("Guess orbital not found: " << orbname.str());
                success &= false;
            }
            if (phi_i.hasReal()) {
                Phi[i].alloc(NUMBER::Real);
                // Refine to get accurate function values
                mrcpp::refine_grid(phi_i.real(), 1);
                mrcpp::project(prec, Phi[i].real(), phi_i.real());
            }
            if (phi_i.hasImag()) {
                Phi[i].alloc(NUMBER::Imag);
                // Refine to get accurate function values
                mrcpp::refine_grid(phi_i.imag(), 1);
                mrcpp::project(prec, Phi[i].imag(), phi_i.imag());
            }
            std::stringstream o_txt;
            o_txt << std::setw(w1 - 1) << i;
            o_txt << std::setw(w3) << print_utils::dbl_to_str(Phi[i].norm(), pprec, true);
            print_utils::qmfunction(1, o_txt.str(), Phi[i], t_i);
        }
    }
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    mrcpp::print::footer(1, t_tot, 2);
    return success;
}

} // namespace mrchem
