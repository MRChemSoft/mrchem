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

#include "mrchem.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/IdentityOperator.h"

using namespace mrchem;
using namespace orbital;

namespace identity_operator_tests {

auto f = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-1.0 * R * R);
};

auto g = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-2.0 * R * R);
};

TEST_CASE("IdentityOperator", "[identity_operator]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("apply") {
        Orbital phi(SPIN::Paired);
        mrcpp::cplxfunc::project(phi, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(phi, g, NUMBER::Imag, prec);

        IdentityOperator I;
        I.setup(prec);
        Orbital Iphi = I(phi);
        I.clear();

        REQUIRE(Iphi.integrate().real() == Catch::Approx(phi.integrate().real()));
        REQUIRE(Iphi.integrate().imag() == Catch::Approx(phi.integrate().imag()));
    }

    SECTION("vector apply") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.distribute();

        if (mrcpp::mpi::my_orb(Phi[0])) mrcpp::cplxfunc::project(Phi[0], f, NUMBER::Real, prec);
        if (mrcpp::mpi::my_orb(Phi[1])) mrcpp::cplxfunc::project(Phi[1], g, NUMBER::Imag, prec);
        normalize(Phi);

        IdentityOperator I;
        I.setup(prec);
        SECTION("O(Phi)") {
            OrbitalVector IPhi = I(Phi);
            ComplexVector ints_a = orbital::get_integrals(Phi);
            ComplexVector ints_b = orbital::get_integrals(IPhi);
            REQUIRE(ints_a.real()(0) == Catch::Approx(ints_b.real()(0)));
            REQUIRE(ints_a.real()(1) == Catch::Approx(ints_b.real()(1)));
            REQUIRE(ints_a.imag()(0) == Catch::Approx(ints_b.imag()(0)));
            REQUIRE(ints_a.imag()(1) == Catch::Approx(ints_b.imag()(1)));
        }
        SECTION("trace") {
            double nEl = get_electron_number(Phi);
            ComplexDouble tr = I.trace(Phi);
            REQUIRE(tr.real() == Catch::Approx(nEl));
            REQUIRE(std::abs(tr.imag()) < thrs);
        }
        I.clear();
    }

    SECTION("expectation value") {
        Orbital phi(SPIN::Paired);
        mrcpp::cplxfunc::project(phi, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(phi, g, NUMBER::Imag, prec);

        IdentityOperator I;
        I.setup(prec);
        ComplexDouble S = I(phi, phi);
        I.clear();

        REQUIRE(S.real() == Catch::Approx(phi.squaredNorm()));
        REQUIRE(S.imag() < thrs);
    }

    SECTION("expectation matrix") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.distribute();

        if (mrcpp::mpi::my_orb(Phi[0])) mrcpp::cplxfunc::project(Phi[0], f, NUMBER::Imag, prec);
        if (mrcpp::mpi::my_orb(Phi[1])) mrcpp::cplxfunc::project(Phi[1], g, NUMBER::Imag, prec);

        IdentityOperator I;
        I.setup(prec);
        ComplexMatrix S = I(Phi, Phi);
        I.clear();

        DoubleMatrix sq_norms = orbital::get_squared_norms(Phi);
        REQUIRE(std::abs(S(0, 0)) == Catch::Approx(sq_norms(0)));
        REQUIRE(std::abs(S(1, 1)) == Catch::Approx(sq_norms(1)));
        REQUIRE(S(0, 1).real() == Catch::Approx(S(1, 0).real()));
        REQUIRE(S(0, 1).imag() == Catch::Approx(-S(1, 0).imag()));
    }
}

} // namespace identity_operator_tests
