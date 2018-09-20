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

#include "qmoperators/one_electron/IdentityOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using namespace mrchem;
using namespace orbital;

namespace identity_operator_tests {

auto f = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-1.0*R*R);
};

auto g = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-2.0*R*R);
};

TEST_CASE("IdentityOperator", "[identity_operator]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("apply") {
        Orbital phi(SPIN::Paired);
        phi.alloc();
        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        IdentityOperator I;
        I.setup(prec);

        Orbital Iphi = I(phi);
        REQUIRE( Iphi.real().integrate() == Approx(phi.real().integrate()) );
        REQUIRE( Iphi.imag().integrate() == Approx(phi.imag().integrate()) );
        Iphi.free();

        I.clear();
        phi.free();
    }

    SECTION("vector apply") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        Phi[0].alloc(NUMBER::Real);
        Phi[1].alloc(NUMBER::Real);
        mrcpp::project(prec, Phi[0].real(), f);
        mrcpp::project(prec, Phi[1].real(), g);
        normalize(Phi);

        IdentityOperator I;
        I.setup(prec);
        SECTION("O(Phi)") {
            OrbitalVector IPhi = I(Phi);
            REQUIRE( IPhi[0].real().integrate() == Approx(Phi[0].real().integrate()) );
            REQUIRE( IPhi[1].real().integrate() == Approx(Phi[1].real().integrate()) );
            free(IPhi);
        }
        SECTION("trace") {
            double nEl = get_electron_number(Phi);
            ComplexDouble tr = I.trace(Phi);
            REQUIRE( tr.real() == Approx(nEl) );
            REQUIRE( std::abs(tr.imag()) < thrs );
        }
        I.clear();
        free(Phi);
    }

    SECTION("expectation value") {
        Orbital phi(SPIN::Paired);
        phi.alloc();
        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        IdentityOperator I;
        I.setup(prec);

        ComplexDouble S = I(phi, phi);
        REQUIRE( S.real() == Approx(phi.squaredNorm()) );
        REQUIRE( S.imag() < thrs );

        I.clear();
        phi.free();
    }

    SECTION("expectation matrix") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        Phi[0].alloc(NUMBER::Imag);
        Phi[1].alloc(NUMBER::Imag);
        mrcpp::project(prec, Phi[0].imag(), f);
        mrcpp::project(prec, Phi[1].imag(), g);

        IdentityOperator I;
        I.setup(prec);

        ComplexMatrix S = I(Phi, Phi);
        REQUIRE( std::abs(S(0,0)) == Approx(Phi[0].squaredNorm()) );
        REQUIRE( std::abs(S(1,1)) == Approx(Phi[1].squaredNorm()) );
        REQUIRE( S(0,1).real() == Approx( S(1,0).real()) );
        REQUIRE( S(0,1).imag() == Approx(-S(1,0).imag()) );

        I.clear();
        free(Phi);
    }
}

} //namespace identity_operator_tests
