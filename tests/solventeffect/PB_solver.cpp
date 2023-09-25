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

#include "catch.hpp"

#include <MRCPP/MWOperators>

#include "mrchem.h"

#include "chemistry/Nucleus.h"
#include "chemistry/PeriodicTable.h"
#include "environment/Cavity.h"
#include "environment/DHScreening.h"
#include "environment/GPESolver.h"
#include "environment/LPBESolver.h"
#include "environment/PBESolver.h"
#include "environment/Permittivity.h"
#include "qmfunctions/qmfunction_fwd.h"
#include "qmoperators/two_electron/ReactionOperator.h"

using namespace mrchem;

namespace PB_solver {

TEST_CASE("Poisson Boltzmann equation solver standard", "[pb_standard]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-8;

    auto acc_pot = true;
    auto dyn_thrs = false;
    auto kain = 7;
    auto max_iter = 100;
    auto eps_in = 1.0;
    auto eps_out = 78.54;
    auto kappa_out = 0.054995;
    auto slope = 0.1;

    auto R = std::vector<double>({3.7794522509156563});
    auto sph_coords = std::vector<mrcpp::Coord<3>>({{0.0, 0.0, 0.0}});
    // initialize spherical cavity
    auto sphere = Cavity(sph_coords, R, slope);
    // initialize dielectric function
    auto dielectric_func = Permittivity(sphere, eps_in, eps_out, "exponential");
    // initialize DHScreening object
    auto kappa_sq = DHScreening(sphere, kappa_out, "continuous");

    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);

    PeriodicTable PT;
    SECTION("case 0: one positive charge in center of sphere (born model)", "[pb_standard][case_0]") {
        // initialize Nuclei  in center of sphere
        auto q_coords = std::vector<mrcpp::Coord<3>>({{0.0, 0.0, 0.0}});
        Nucleus Q(PT.getElement(0), q_coords[0]);
        Nuclei molecule;
        molecule.push_back(Q);

        // initialize orbital vector
        auto Phi_p = std::make_shared<OrbitalVector>();
        auto &Phi = *Phi_p;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.distribute();

        auto scrf_p = std::make_unique<PBESolver>(dielectric_func, kappa_sq, molecule, P_p, D_p, prec / 100, kain, max_iter, acc_pot, dyn_thrs, "nuclear");
        auto Reo = std::make_shared<ReactionOperator>(std::move(scrf_p), Phi_p);
        Reo->setup(prec * 10);

        double E_s = Reo->getTotalEnergy() / 2.0;
        Reo->clear();
        REQUIRE(E_s == Approx(-0.1373074208).epsilon(thrs));
    }
}

TEST_CASE("Poisson Boltzmann equation solver linearized", "[pb_linearized]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-8;

    auto acc_pot = true;
    auto dyn_thrs = false;
    auto kain = 5;
    auto max_iter = 200;
    auto eps_in = 1.0;
    auto eps_out = 78.54;
    auto kappa_out = 0.054995;
    auto slope = 0.1;

    auto R = std::vector<double>({3.7794522509156563});
    auto sph_coords = std::vector<mrcpp::Coord<3>>({{0.0, 0.0, 0.0}});
    // initialize spherical cavity
    auto sphere = Cavity(sph_coords, R, slope);
    // initialize dielectric function
    auto dielectric_func = Permittivity(sphere, eps_in, eps_out, "exponential");
    // initialize DHScreening object
    auto kappa_sq = DHScreening(sphere, kappa_out, "continuous");

    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);

    PeriodicTable PT;
    SECTION("case 0: one positive charge in center of sphere (born model)", "[pb_linearized][case_0]") {

        // initialize Nuclei  in center of sphere
        auto q_coords = std::vector<mrcpp::Coord<3>>({{0.0, 0.0, 0.0}});
        Nucleus Q(PT.getElement(0), q_coords[0]);
        Nuclei molecule;
        molecule.push_back(Q);

        // initialize orbital vector
        auto Phi_p = std::make_shared<OrbitalVector>();
        auto &Phi = *Phi_p;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.distribute();

        auto scrf_p = std::make_unique<LPBESolver>(dielectric_func, kappa_sq, molecule, P_p, D_p, prec / 100, kain, max_iter, acc_pot, dyn_thrs, "nuclear");

        auto Reo = std::make_shared<ReactionOperator>(std::move(scrf_p), Phi_p);
        Reo->setup(prec * 10);

        double E_s = Reo->getTotalEnergy() / 2.0;
        Reo->clear();
        REQUIRE(E_s == Approx(-0.1329646842).epsilon(thrs));
    }
}

} // namespace PB_solver
