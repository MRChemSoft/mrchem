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

using namespace mrchem;

namespace orbital_tests {

auto f = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-1.0 * R * R);
};

auto g = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-2.0 * R * R);
};

TEST_CASE("Orbital", "[orbital]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("copy orbital") {
        Orbital phi_1(SPIN::Paired);
        mrcpp::cplxfunc::project(phi_1, f, NUMBER::Real, prec);

        SECTION("copy constructor") {
            Orbital phi_2(phi_1);
            REQUIRE(phi_2.occ() == phi_1.occ());
            REQUIRE(phi_2.spin() == phi_1.spin());
            REQUIRE(phi_2.norm() == phi_1.norm());
            REQUIRE(&phi_2.real() == &phi_1.real());
            REQUIRE(&phi_2.imag() == &phi_1.imag());
        }

        SECTION("default constructor plus assignment") {
            Orbital phi_2;
            phi_2 = phi_1;
            REQUIRE(phi_2.occ() == phi_1.occ());
            REQUIRE(phi_2.spin() == phi_1.spin());
            REQUIRE(phi_2.norm() == phi_1.norm());
            REQUIRE(&phi_2.real() == &phi_1.real());
            REQUIRE(&phi_2.imag() == &phi_1.imag());
        }

        SECTION("assigment constructor") {
            Orbital phi_2 = phi_1;
            REQUIRE(phi_2.occ() == phi_1.occ());
            REQUIRE(phi_2.spin() == phi_1.spin());
            REQUIRE(phi_2.norm() == phi_1.norm());
            REQUIRE(&phi_2.real() == &phi_1.real());
            REQUIRE(&phi_2.imag() == &phi_1.imag());
        }

        SECTION("deep copy") {
            Orbital phi_2(SPIN::Alpha);
            mrcpp::cplxfunc::deep_copy(phi_2, phi_1);
            REQUIRE(phi_2.occ() != phi_1.occ());
            REQUIRE(phi_2.spin() != phi_1.spin());
            REQUIRE(phi_2.norm() == phi_1.norm());
            REQUIRE(&phi_2.real() != &phi_1.real());
            REQUIRE(not(phi_2.hasImag()));
        }

        SECTION("parameter copy") {
            Orbital phi_2 = phi_1.paramCopy();
            REQUIRE(phi_2.occ() == phi_1.occ());
            REQUIRE(phi_2.spin() == phi_1.spin());
            REQUIRE(phi_2.norm() < 1.0);
            REQUIRE(not(phi_2.hasReal()));
            REQUIRE(not(phi_2.hasImag()));
        }
    }

    SECTION("normalize") {
        Orbital phi(SPIN::Paired);
        REQUIRE(phi.norm() == Catch::Approx(-1.0));

        mrcpp::cplxfunc::project(phi, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(phi, g, NUMBER::Imag, prec);
        REQUIRE(phi.norm() > 1.0);

        orbital::normalize(phi);
        REQUIRE(phi.norm() == Catch::Approx(1.0));
    }

    SECTION("orthogonalize") {
        Orbital phi_1(SPIN::Alpha);
        mrcpp::cplxfunc::project(phi_1, f, NUMBER::Real, prec);

        WHEN("orbitals have different spins") {
            Orbital phi_2(SPIN::Beta);
            mrcpp::cplxfunc::project(phi_2, g, NUMBER::Imag, prec);

            THEN("their overlap is zero") {
                ComplexDouble S = orbital::dot(phi_1, phi_2);
                REQUIRE(std::abs(S.real()) < thrs);
                REQUIRE(std::abs(S.imag()) < thrs);
            }
        }

        WHEN("orbitals have the same spin") {
            Orbital phi_2(SPIN::Alpha);
            mrcpp::cplxfunc::project(phi_2, g, NUMBER::Imag, prec);

            THEN("their overlap is non-zero") {
                ComplexDouble S1 = orbital::dot(phi_1, phi_2);
                REQUIRE(std::abs(S1.real()) < thrs);
                REQUIRE(std::abs(S1.imag()) > thrs);

                AND_THEN("<phi_1|phi_2^dag> = <phi_1|phi_2>*") {
                    ComplexDouble S2 = orbital::dot(phi_1, phi_2.dagger());
                    REQUIRE(S2.real() == Catch::Approx(S1.real()));
                    REQUIRE(S2.imag() == Catch::Approx(-S1.imag()));
                }
            }

            AND_THEN("they are orthogonalized") {
                orbital::orthogonalize(prec, phi_2, phi_1);

                THEN("their overlap is zero") {
                    ComplexDouble S3 = orbital::dot(phi_1, phi_2);
                    REQUIRE(std::abs(S3.real()) < thrs);
                    REQUIRE(std::abs(S3.imag()) < thrs);
                }
            }
        }
    }
}

} // namespace orbital_tests
