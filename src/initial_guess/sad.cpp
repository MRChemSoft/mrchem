/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <MRCPP/utils/details.h>

#include "core.h"
#include "gto.h"
#include "parallel.h"
#include "sad.h"

#include "chemistry/Nucleus.h"

#include "utils/print_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

#include "mrdft/Factory.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace sad {

void project_atomic_densities(double prec, Density &rho_tot, const Nuclei &nucs);

} // namespace sad
} // namespace initial_guess

bool initial_guess::sad::setup(OrbitalVector &Phi, double prec, const Nuclei &nucs, int zeta) {
    if (Phi.size() == 0) return false;

    auto restricted = (orbital::size_singly(Phi)) ? false : true;
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Compute initial orbitals");
    print_utils::text(0, "Method      ", "Diagonalize SAD Hamiltonian");
    print_utils::text(0, "Precision   ", print_utils::dbl_to_str(prec, 5, true));
    print_utils::text(0, "Restricted  ", (restricted) ? "True" : "False");
    print_utils::text(0, "Functional  ", "LDA (SVWN5)");
    print_utils::text(0, "AO basis    ", "Hydrogenic orbitals");
    print_utils::text(0, "Zeta quality", std::to_string(zeta));
    mrcpp::print::separator(0, '~', 2);

    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "SAD Initial Guess");

    // Make Fock operator contributions
    t_lap.start();
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);

    mrdft::Factory xc_factory(*MRA);
    xc_factory.setSpin(false);
    xc_factory.setFunctional("SLATERX", 1.0);
    xc_factory.setFunctional("VWN5C", 1.0);
    auto mrdft_p = xc_factory.build();

    MomentumOperator p(D_p);
    NuclearOperator V_nuc(nucs, prec);
    CoulombOperator J(P_p);
    XCOperator XC(mrdft_p);
    RankZeroOperator V = V_nuc + J + XC;
    if (plevel == 1) mrcpp::print::time(1, "Projecting nuclear potential", t_lap);

    // Compute Coulomb density
    t_lap.start();
    Density &rho_j = J.getDensity();
    initial_guess::sad::project_atomic_densities(prec, rho_j, nucs);

    // Compute XC density
    Density &rho_xc = XC.getDensity(DensityType::Total);
    qmfunction::deep_copy(rho_xc, rho_j);
    if (plevel == 1) mrcpp::print::time(1, "Projecting GTO density", t_lap);

    // Project AO basis of hydrogen functions
    t_lap.start();
    OrbitalVector Psi;
    initial_guess::core::project_ao(Psi, prec, nucs, zeta);
    if (plevel == 1) mrcpp::print::time(1, "Projecting Hydrogen AOs", t_lap);

    mrcpp::print::header(2, "Building Fock operator");
    t_lap.start();
    p.setup(prec);
    V.setup(prec);
    mrcpp::print::footer(2, t_lap, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_lap);

    // Compute Fock matrix
    mrcpp::print::header(2, "Diagonalizing Fock matrix");
    ComplexMatrix U = initial_guess::core::diagonalize(Psi, p, V);

    // Rotate orbitals and fill electrons by Aufbau
    t_lap.start();
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);
    initial_guess::core::rotate_orbitals(Phi, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_a, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_b, prec, U, Psi);
    for (auto &phi_a : Phi_a) Phi.push_back(phi_a);
    for (auto &phi_b : Phi_b) Phi.push_back(phi_b);
    p.clear();
    V.clear();

    mrcpp::print::footer(2, t_lap, 2);
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
    return true;
}

void initial_guess::sad::project_atomic_densities(double prec, Density &rho_tot, const Nuclei &nucs) {
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 8;
    auto w3 = w0 / 3;
    auto w4 = w0 - (w1 + w2 + 2 * w3);

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::setw(w4) << " ";
    o_head << std::setw(w3) << "Nuclear charge";
    o_head << std::setw(w3) << "Electron charge";

    mrcpp::print::header(2, "Projecting GTO density");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    auto crop_prec = (mpi::numerically_exact) ? -1.0 : prec;
    std::string sad_path;
    for (auto n : {sad_basis_source_dir(), sad_basis_install_dir()}) {
        auto trimmed = print_utils::rtrim_copy(n);
        if (mrcpp::details::directory_exists(trimmed)) {
            sad_path = trimmed;
            break;
        }
    }

    Timer t_tot;
    Density rho_loc(false);
    rho_loc.alloc(NUMBER::Real);
    rho_loc.real().setZero();

    auto tot_nuc = 0.0;
    auto tot_rho = 0.0;

    Timer t_loc;
    for (int k = 0; k < nucs.size(); k++) {
        if (mpi::orb_rank != k % mpi::orb_size) continue;

        const std::string &sym = nucs[k].getElement().getSymbol();
        std::stringstream o_bas, o_dens;
        o_bas << sad_path << "/" << sym << ".bas";
        o_dens << sad_path << "/" << sym << ".dens";

        Density rho_k = initial_guess::gto::project_density(prec, nucs[k], o_bas.str(), o_dens.str());
        rho_loc.add(1.0, rho_k);
        rho_loc.crop(crop_prec);

        auto nuc_charge = nucs[k].getCharge();
        auto rho_charge = rho_k.integrate().real();
        tot_nuc += nuc_charge;
        tot_rho += rho_charge;

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << sym;
        o_row << std::setw(w4) << " ";
        o_row << std::setw(w3) << print_utils::dbl_to_str(nuc_charge, 2 * pprec, false);
        o_row << std::setw(w3) << print_utils::dbl_to_str(rho_charge, 2 * pprec, false);
        println(2, o_row.str());
    }
    t_loc.stop();

    Timer t_com;
    density::allreduce_density(prec, rho_tot, rho_loc);
    t_com.stop();

    std::stringstream o_row;
    o_row << " Total charge";
    o_row << std::string(w1 + w2 + w4 - 13, ' ');
    o_row << std::setw(w3) << print_utils::dbl_to_str(tot_nuc, 2 * pprec, false);
    o_row << std::setw(w3) << print_utils::dbl_to_str(tot_rho, 2 * pprec, false);

    mrcpp::print::separator(2, '-');
    println(2, o_row.str());
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local density", rho_loc, t_loc);
    print_utils::qmfunction(2, "Allreduce density", rho_tot, t_com);
    mrcpp::print::footer(2, t_tot, 2);
}

} // namespace mrchem
