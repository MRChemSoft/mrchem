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

#include "catch2/catch_all.hpp"

#include "analyticfunctions/HydrogenFunction.h"
#include "mrchem.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"

using namespace mrchem;

namespace density_tests {

TEST_CASE("Density", "[density]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("calc density") {
        OrbitalVector Phi;
        for (int i = 0; i < 5; i++) Phi.push_back(Orbital(SPIN::Alpha));
        for (int i = 0; i < 2; i++) Phi.push_back(Orbital(SPIN::Beta));
        Phi.distribute();

        HydrogenFunction s1(1, 0, 0);
        HydrogenFunction s2(2, 0, 0);
        HydrogenFunction px(2, 1, 0);
        HydrogenFunction py(2, 1, 1);
        HydrogenFunction pz(2, 1, 2);
        if (mrcpp::mpi::my_orb(Phi[0])) mrcpp::cplxfunc::project(Phi[0], s1, NUMBER::Real, prec);
        if (mrcpp::mpi::my_orb(Phi[1])) mrcpp::cplxfunc::project(Phi[1], s2, NUMBER::Real, prec);
        if (mrcpp::mpi::my_orb(Phi[2])) mrcpp::cplxfunc::project(Phi[2], px, NUMBER::Imag, prec);
        if (mrcpp::mpi::my_orb(Phi[3])) mrcpp::cplxfunc::project(Phi[3], py, NUMBER::Imag, prec);
        if (mrcpp::mpi::my_orb(Phi[4])) mrcpp::cplxfunc::project(Phi[4], pz, NUMBER::Imag, prec);
        if (mrcpp::mpi::my_orb(Phi[5])) mrcpp::cplxfunc::project(Phi[5], s1, NUMBER::Real, prec);
        if (mrcpp::mpi::my_orb(Phi[6])) mrcpp::cplxfunc::project(Phi[6], s1, NUMBER::Real, prec);

        SECTION("non-shared memory total/spin density") {
            Density rho_t(false);
            Density rho_s(false);

            density::compute(prec, rho_t, Phi, DensityType::Total);
            density::compute(prec, rho_s, Phi, DensityType::Spin);

            REQUIRE(rho_t.integrate().real() == Catch::Approx(7.0));
            REQUIRE(rho_s.integrate().real() == Catch::Approx(3.0));
        }

        SECTION("shared memory alpha/beta density") {
            Density rho_a(true);
            Density rho_b(true);

            density::compute(prec, rho_a, Phi, DensityType::Alpha);
            density::compute(prec, rho_b, Phi, DensityType::Beta);

            REQUIRE(rho_a.integrate().real() == Catch::Approx(5.0));
            REQUIRE(rho_b.integrate().real() == Catch::Approx(2.0));
        }
    }
}

} // namespace density_tests
