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
#include "utils/math_utils.h"

using namespace mrchem;
using namespace orbital;

namespace orbital_vector_tests {

auto f1 = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-1.0 * R * R);
};

auto f2 = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-2.0 * R * R);
};

auto f3 = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-3.0 * R * R);
};

auto f4 = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-4.0 * R * R);
};

auto f5 = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-5.0 * R * R);
};

auto f6 = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-6.0 * R * R);
};

TEST_CASE("OrbitalVector", "[orbital_vector]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("push_back") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Beta));
        Phi.push_back(Orbital(SPIN::Alpha));
        Phi.push_back(Orbital(SPIN::Beta));

        REQUIRE(Phi.size() == 4);
        REQUIRE(get_electron_number(Phi, SPIN::Paired) == 5);
        REQUIRE(get_electron_number(Phi, SPIN::Alpha) == 2);
        REQUIRE(get_electron_number(Phi, SPIN::Beta) == 3);
        REQUIRE(size_empty(Phi) == 0);
        REQUIRE(size_paired(Phi) == 1);
        REQUIRE(size_alpha(Phi) == 1);
        REQUIRE(size_beta(Phi) == 2);
        REQUIRE(get_multiplicity(Phi) == 2);

        Phi.clear();
        REQUIRE(Phi.size() == 0);
    }

    SECTION("adjoin/disjoin vectors") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Alpha));
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Beta));
        Phi.push_back(Orbital(SPIN::Beta));
        Phi.distribute();

        OrbitalVector Phi_p = disjoin(Phi, SPIN::Paired);
        OrbitalVector Phi_a = disjoin(Phi, SPIN::Alpha);
        OrbitalVector Phi_b = disjoin(Phi, SPIN::Beta);

        REQUIRE(Phi.size() == 0);
        REQUIRE(Phi_p.size() == 2);
        REQUIRE(Phi_a.size() == 1);
        REQUIRE(Phi_b.size() == 2);

        Phi = adjoin(Phi, Phi_p);
        Phi = adjoin(Phi, Phi_a);
        Phi = adjoin(Phi, Phi_b);

        REQUIRE(Phi.size() == 5);
        REQUIRE(Phi_p.size() == 0);
        REQUIRE(Phi_a.size() == 0);
        REQUIRE(Phi_b.size() == 0);

        REQUIRE(Phi[0].spin() == SPIN::Paired);
        REQUIRE(Phi[1].spin() == SPIN::Paired);
        REQUIRE(Phi[2].spin() == SPIN::Alpha);
        REQUIRE(Phi[3].spin() == SPIN::Beta);
        REQUIRE(Phi[4].spin() == SPIN::Beta);
    }

    SECTION("copy vectors") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Alpha));
        Phi.distribute();
        if (mrcpp::mpi::my_orb(Phi[0])) mrcpp::cplxfunc::project(Phi[0], f1, NUMBER::Real, prec);
        if (mrcpp::mpi::my_orb(Phi[1])) mrcpp::cplxfunc::project(Phi[1], f2, NUMBER::Imag, prec);
        normalize(Phi);

        SECTION("copy constructor") {
            OrbitalVector Psi(Phi);
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Catch::Approx(1.0));
            REQUIRE(norms[1] == Catch::Approx(1.0));
        }

        SECTION("default constructor plus assignment") {
            OrbitalVector Psi;
            Psi = Phi;
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Catch::Approx(1.0));
            REQUIRE(norms[1] == Catch::Approx(1.0));
        }

        SECTION("default constructor plus deep copy") {
            OrbitalVector Psi;
            Psi = deep_copy(Phi);
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Catch::Approx(1.0));
            REQUIRE(norms[1] == Catch::Approx(1.0));
        }

        SECTION("assigment constructor") {
            OrbitalVector Psi = Phi;
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Catch::Approx(1.0));
            REQUIRE(norms[1] == Catch::Approx(1.0));
        }

        SECTION("parameter copy") {
            OrbitalVector Psi;
            Psi = param_copy(Phi);
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] < 0.0);
            REQUIRE(norms[1] < 0.0);
        }
    }

    SECTION("normalization") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Alpha));
        Phi.distribute();

        DoubleVector norms1 = get_norms(Phi);
        REQUIRE(norms1[0] == Catch::Approx(-1.0));
        REQUIRE(norms1[1] == Catch::Approx(-1.0));

        if (mrcpp::mpi::my_orb(Phi[0])) mrcpp::cplxfunc::project(Phi[0], f1, NUMBER::Real, prec);
        if (mrcpp::mpi::my_orb(Phi[1])) mrcpp::cplxfunc::project(Phi[1], f2, NUMBER::Real, prec);

        DoubleVector norms2 = get_norms(Phi);
        REQUIRE(norms2[0] > 0.0);
        REQUIRE(norms2[1] > 0.0);

        normalize(Phi);

        DoubleVector norms3 = get_norms(Phi);
        REQUIRE(norms3[0] == Catch::Approx(1.0));
        REQUIRE(norms3[1] == Catch::Approx(1.0));

        SECTION("norm overlap") {
            ComplexMatrix S = orbital::calc_overlap_matrix(Phi);
            DoubleMatrix Snorm = orbital::calc_norm_overlap_matrix(Phi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(Snorm(i, j) == Catch::Approx(1.0));
                    if (i != j) REQUIRE(Snorm(i, j) >= Catch::Approx(std::abs(S(i, j))));
                }
            }
        }
    }

    SECTION("orthogonalization") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Beta));
        Phi.push_back(Orbital(SPIN::Alpha));
        Phi.push_back(Orbital(SPIN::Alpha));
        Phi.push_back(Orbital(SPIN::Beta));
        Phi.distribute();

        if (true or mrcpp::mpi::my_orb(Phi[0])) mrcpp::cplxfunc::project(Phi[0], f1, NUMBER::Real, prec);
        if (true or mrcpp::mpi::my_orb(Phi[1])) mrcpp::cplxfunc::project(Phi[1], f2, NUMBER::Real, prec);
        if (true or mrcpp::mpi::my_orb(Phi[2])) mrcpp::cplxfunc::project(Phi[2], f3, NUMBER::Real, prec);
        if (true or mrcpp::mpi::my_orb(Phi[3])) mrcpp::cplxfunc::project(Phi[3], f4, NUMBER::Real, prec);

        // Complex phase rotation
        for (int n = 0; n < Phi.size(); n++) {
            double theta = (n + 1.0) * mrcpp::pi / 7.0;
            ComplexDouble phase(std::cos(theta), std::sin(theta));
            Phi[n].rescale(phase);
        }

        SECTION("Lowdin orthonormalize") {
            ComplexMatrix M = ComplexMatrix::Zero(Phi.size(), Phi.size());
            orthonormalize(-1.0, Phi, M);

            ComplexMatrix S = calc_overlap_matrix(Phi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(std::abs(S(i, j)) == Catch::Approx(1.0));
                    if (i != j) REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
        }

        SECTION("diagonalize overlap") {
            auto S = calc_overlap_matrix(Phi);
            diagonalize(-1.0, Phi, S);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(std::abs(S(i, j)) == Catch::Approx(1.0));
                    if (i != j) REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
        }

        SECTION("Gram-Schmidt orthonormalize") {
            orthogonalize(prec, Phi);
            normalize(Phi);

            ComplexMatrix S = calc_overlap_matrix(Phi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(std::abs(S(i, j)) == Catch::Approx(1.0));
                    if (i != j) REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
        }

        SECTION("vector orthogonalize") {
            OrbitalVector Psi;
            Psi.push_back(Orbital(SPIN::Alpha));
            Psi.push_back(Orbital(SPIN::Beta));
            Psi.distribute();

            if (mrcpp::mpi::my_orb(Psi[0])) mrcpp::cplxfunc::project(Psi[0], f5, NUMBER::Real, prec);
            if (mrcpp::mpi::my_orb(Psi[1])) mrcpp::cplxfunc::project(Psi[1], f6, NUMBER::Real, prec);

            orthogonalize(prec, Phi);
            orthogonalize(prec, Psi, Phi);

            ComplexMatrix S = orbital::calc_overlap_matrix(Psi, Phi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) { REQUIRE(std::abs(S(i, j)) < thrs); }
            }
        }
    }

    SECTION("orbital transformations") {
        OrbitalVector Phi;
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.push_back(Orbital(SPIN::Paired));
        Phi.distribute();

        if (mrcpp::mpi::my_orb(Phi[0])) {
            mrcpp::cplxfunc::project(Phi[0], f1, NUMBER::Real, prec);
            mrcpp::cplxfunc::project(Phi[0], f2, NUMBER::Imag, prec);
        }

        if (mrcpp::mpi::my_orb(Phi[1])) {
            mrcpp::cplxfunc::project(Phi[1], f3, NUMBER::Real, prec);
            mrcpp::cplxfunc::project(Phi[1], f4, NUMBER::Imag, prec);
        }

        orthogonalize(prec, Phi);
        normalize(Phi);

        SECTION("unitary transformation") {
            double theta = 0.5;
            ComplexMatrix U(Phi.size(), Phi.size());
            U(0, 0) = std::cos(theta);
            U(0, 1) = -std::sin(theta);
            U(1, 0) = std::sin(theta);
            U(1, 1) = std::cos(theta);

            OrbitalVector Psi = rotate(Phi, U, prec);
            ComplexMatrix S = orbital::calc_overlap_matrix(Psi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(std::abs(S(i, j)) == Catch::Approx(1.0));
                    if (i != j) REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
        }

        SECTION("vector addition") {
            // Complex phase rotation
            double theta = 0.6;
            ComplexDouble a(std::cos(theta), 0.0), b(0.0, std::sin(theta));
            OrbitalVector Psi = add(a, Phi, b, Phi, prec);

            ComplexMatrix S = orbital::calc_overlap_matrix(Psi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(std::abs(S(i, j)) == Catch::Approx(1.0));
                    if (i != j) REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
        }
    }
}

} // namespace orbital_vector_tests
