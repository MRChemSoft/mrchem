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

#include "mrchem.h"
#include "parallel.h"

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
        phi.alloc(NUMBER::Total);
        mrcpp::project<3>(prec, phi.real(), f);
        mrcpp::project<3>(prec, phi.imag(), g);

        IdentityOperator I;
        I.setup(prec);

        Orbital Iphi = I(phi);
        REQUIRE(Iphi.real().integrate() == Approx(phi.real().integrate()));
        REQUIRE(Iphi.imag().integrate() == Approx(phi.imag().integrate()));
        Iphi.free();

        I.clear();
        phi.free();
    }

    SECTION("vector apply") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        mpi::distribute(Phi);

        if (mpi::my_orb(Phi[0])) Phi[0].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[1])) Phi[1].alloc(NUMBER::Imag);
        if (mpi::my_orb(Phi[0])) mrcpp::project<3>(prec, Phi[0].real(), f);
        if (mpi::my_orb(Phi[1])) mrcpp::project<3>(prec, Phi[1].imag(), g);
        normalize(Phi);

        IdentityOperator I;
        I.setup(prec);
        SECTION("O(Phi)") {
            OrbitalVector IPhi = I(Phi);
            ComplexVector ints_a = orbital::get_integrals(Phi);
            ComplexVector ints_b = orbital::get_integrals(IPhi);
            REQUIRE(ints_a.real()(0) == Approx(ints_b.real()(0)));
            REQUIRE(ints_a.real()(1) == Approx(ints_b.real()(1)));
            REQUIRE(ints_a.imag()(0) == Approx(ints_b.imag()(0)));
            REQUIRE(ints_a.imag()(1) == Approx(ints_b.imag()(1)));
            free(IPhi);
        }
        SECTION("trace") {
            double nEl = get_electron_number(Phi);
            ComplexDouble tr = I.trace(Phi);
            REQUIRE(tr.real() == Approx(nEl));
            REQUIRE(std::abs(tr.imag()) < thrs);
        }
        I.clear();
        free(Phi);
    }

    SECTION("expectation value") {
        Orbital phi(SPIN::Paired);
        phi.alloc(NUMBER::Total);
        mrcpp::project<3>(prec, phi.real(), f);
        mrcpp::project<3>(prec, phi.imag(), g);

        IdentityOperator I;
        I.setup(prec);

        ComplexDouble S = I(phi, phi);
        REQUIRE(S.real() == Approx(phi.squaredNorm()));
        REQUIRE(S.imag() < thrs);

        I.clear();
        phi.free();
    }

    SECTION("expectation matrix") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        mpi::distribute(Phi);

        if (mpi::my_orb(Phi[0])) Phi[0].alloc(NUMBER::Imag);
        if (mpi::my_orb(Phi[1])) Phi[1].alloc(NUMBER::Imag);
        if (mpi::my_orb(Phi[0])) mrcpp::project<3>(prec, Phi[0].imag(), f);
        if (mpi::my_orb(Phi[1])) mrcpp::project<3>(prec, Phi[1].imag(), g);

        IdentityOperator I;
        I.setup(prec);

        ComplexMatrix S = I(Phi, Phi);
        DoubleMatrix sq_norms = orbital::get_squared_norms(Phi);

        REQUIRE(std::abs(S(0, 0)) == Approx(sq_norms(0)));
        REQUIRE(std::abs(S(1, 1)) == Approx(sq_norms(1)));
        REQUIRE(S(0, 1).real() == Approx(S(1, 0).real()));
        REQUIRE(S(0, 1).imag() == Approx(-S(1, 0).imag()));

        I.clear();
        free(Phi);
    }
}

} //namespace identity_operator_tests
