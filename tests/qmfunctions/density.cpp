/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "analyticfunctions/HydrogenFunction.h"
#include "mrchem.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

using namespace mrchem;

namespace density_tests {

TEST_CASE("Density", "[density]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("calc density") {
        Density rho(false);

        SECTION("orbital vector") {
            HydrogenFunction h_1(2, 1, 0);
            HydrogenFunction h_2(2, 1, 1);
            HydrogenFunction h_3(2, 1, 2);

            OrbitalVector Phi;
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Paired);
            Phi.push_back(SPIN::Alpha);
            Phi[0].alloc(NUMBER::Real);
            Phi[1].alloc(NUMBER::Real);
            Phi[2].alloc(NUMBER::Imag);
            mrcpp::project(prec, Phi[0].real(), h_1);
            mrcpp::project(prec, Phi[1].real(), h_2);
            mrcpp::project(prec, Phi[2].imag(), h_3);

            density::compute(prec, rho, Phi, DENSITY::Total);
            REQUIRE(rho.real().integrate() == Approx(4.0));
            orbital::free(Phi);
        }
    }

    SECTION("calc spin density") {
        Density rho_t(false);
        Density rho_s(false);
        Density rho_a(false);
        Density rho_b(false);

        SECTION("orbital vector") {
            HydrogenFunction s1(1, 0, 0);
            HydrogenFunction s2(2, 0, 0);
            HydrogenFunction px(2, 1, 0);
            HydrogenFunction py(2, 1, 1);
            HydrogenFunction pz(2, 1, 2);

            OrbitalVector Phi;
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Beta);
            Phi.push_back(SPIN::Beta);

            Phi[0].alloc(NUMBER::Real);
            Phi[1].alloc(NUMBER::Real);
            Phi[2].alloc(NUMBER::Imag);
            Phi[3].alloc(NUMBER::Imag);
            Phi[4].alloc(NUMBER::Imag);
            Phi[5].alloc(NUMBER::Real);
            Phi[6].alloc(NUMBER::Real);

            mrcpp::project(prec, Phi[0].real(), s1);
            mrcpp::project(prec, Phi[1].real(), s2);
            mrcpp::project(prec, Phi[2].imag(), px);
            mrcpp::project(prec, Phi[3].imag(), py);
            mrcpp::project(prec, Phi[4].imag(), pz);
            mrcpp::project(prec, Phi[5].real(), s1);
            mrcpp::project(prec, Phi[6].real(), s2);

            density::compute(prec, rho_t, Phi, DENSITY::Total);
            density::compute(prec, rho_s, Phi, DENSITY::Spin);
            density::compute(prec, rho_a, Phi, DENSITY::Alpha);
            density::compute(prec, rho_b, Phi, DENSITY::Beta);

            REQUIRE(rho_t.real().integrate() == Approx(7.0));
            REQUIRE(rho_s.real().integrate() == Approx(3.0));
            REQUIRE(rho_a.real().integrate() == Approx(5.0));
            REQUIRE(rho_b.real().integrate() == Approx(2.0));

            orbital::free(Phi);
        }
    }
}

} // namespace density_tests
