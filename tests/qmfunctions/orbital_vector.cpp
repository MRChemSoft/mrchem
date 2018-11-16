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
        Orbital phi_a(SPIN::Alpha);
        Orbital phi_b(SPIN::Beta);

        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(phi_b);
        Phi.push_back(phi_a);
        Phi.push_back(phi_b);

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

    SECTION("alloc") {
        Orbital phi_a(SPIN::Alpha);
        Orbital phi_b(SPIN::Beta);

        phi_a.alloc(NUMBER::Real);
        phi_b.alloc(NUMBER::Real);
        phi_a.real().setZero();
        phi_b.real().setZero();

        OrbitalVector Phi;
        SECTION("clear") {
            Phi.push_back(phi_a);
            Phi.push_back(phi_b);
            Phi.push_back(phi_a);
            Phi.clear();
            phi_a.free();
            phi_b.free();
        }
        SECTION("free") {
            Phi.push_back(phi_a);
            Phi.push_back(phi_b);
            Phi.push_back(phi_a.deepCopy());
            free(Phi);
        }
        REQUIRE(Phi.size() == 0);
    }

    SECTION("adjoin/disjoin vectors") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Alpha);
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Beta);
        Phi.push_back(SPIN::Beta);
        mpi::distribute(Phi);

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
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Alpha);
        mpi::distribute(Phi);

        if (mpi::my_orb(Phi[0])) Phi[0].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[1])) Phi[1].alloc(NUMBER::Imag);
        if (mpi::my_orb(Phi[0])) mrcpp::project<3>(prec, Phi[0].real(), f1);
        if (mpi::my_orb(Phi[1])) mrcpp::project<3>(prec, Phi[1].imag(), f2);
        normalize(Phi);

        SECTION("copy constructor") {
            OrbitalVector Psi(Phi);
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Approx(1.0));
            REQUIRE(norms[1] == Approx(1.0));
            Psi.clear();
        }

        SECTION("default constructor plus assignment") {
            OrbitalVector Psi;
            Psi = Phi;
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Approx(1.0));
            REQUIRE(norms[1] == Approx(1.0));
            Psi.clear();
        }

        SECTION("default constructor plus deep copy") {
            OrbitalVector Psi;
            Psi = deep_copy(Phi);
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Approx(1.0));
            REQUIRE(norms[1] == Approx(1.0));
            free(Psi);
        }

        SECTION("assigment constructor") {
            OrbitalVector Psi = Phi;
            REQUIRE(get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired));
            REQUIRE(get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha));
            REQUIRE(get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta));

            DoubleVector norms = get_norms(Psi);
            REQUIRE(norms[0] == Approx(1.0));
            REQUIRE(norms[1] == Approx(1.0));
            Psi.clear();
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
            Psi.clear();
        }

        free(Phi);
    }

    SECTION("normalization") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Alpha);
        mpi::distribute(Phi);

        DoubleVector norms1 = get_norms(Phi);
        REQUIRE(norms1[0] == Approx(-1.0));
        REQUIRE(norms1[1] == Approx(-1.0));

        if (mpi::my_orb(Phi[0])) Phi[0].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[1])) Phi[1].alloc(NUMBER::Imag);
        if (mpi::my_orb(Phi[0])) mrcpp::project<3>(prec, Phi[0].real(), f1);
        if (mpi::my_orb(Phi[1])) mrcpp::project<3>(prec, Phi[1].imag(), f2);

        DoubleVector norms2 = get_norms(Phi);
        REQUIRE(norms2[0] > 0.0);
        REQUIRE(norms2[1] > 0.0);

        normalize(Phi);

        DoubleVector norms3 = get_norms(Phi);
        REQUIRE(norms3[0] == Approx(1.0));
        REQUIRE(norms3[1] == Approx(1.0));

        free(Phi);
    }

    SECTION("orthogonalization") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Beta);
        Phi.push_back(SPIN::Alpha);
        Phi.push_back(SPIN::Alpha);
        Phi.push_back(SPIN::Beta);
        mpi::distribute(Phi);

        if (mpi::my_orb(Phi[0])) Phi[0].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[1])) Phi[1].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[2])) Phi[2].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[3])) Phi[3].alloc(NUMBER::Real);

        if (mpi::my_orb(Phi[0])) mrcpp::project<3>(prec, Phi[0].real(), f1);
        if (mpi::my_orb(Phi[1])) mrcpp::project<3>(prec, Phi[1].real(), f2);
        if (mpi::my_orb(Phi[2])) mrcpp::project<3>(prec, Phi[2].real(), f3);
        if (mpi::my_orb(Phi[3])) mrcpp::project<3>(prec, Phi[3].real(), f4);

        orthogonalize(Phi);

        SECTION("in place orthonormalize") {
            normalize(Phi);
            ComplexMatrix S = orbital::calc_overlap_matrix(Phi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(std::abs(S(i, j)) == Approx(1.0));
                    if (i != j) REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
        }

        SECTION("vector orthogonalize") {
            OrbitalVector Psi;
            Psi.push_back(SPIN::Alpha);
            Psi.push_back(SPIN::Beta);
            mpi::distribute(Psi);

            if (mpi::my_orb(Psi[0])) Psi[0].alloc(NUMBER::Real);
            if (mpi::my_orb(Psi[1])) Psi[1].alloc(NUMBER::Real);
            if (mpi::my_orb(Psi[0])) mrcpp::project<3>(prec, Psi[0].real(), f5);
            if (mpi::my_orb(Psi[1])) mrcpp::project<3>(prec, Psi[1].real(), f6);

            orthogonalize(Psi, Phi);

            ComplexMatrix S = orbital::calc_overlap_matrix(Psi, Phi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
            free(Psi);
        }
        free(Phi);
    }

    SECTION("addition") {
        OrbitalVector Phi_a;
        Phi_a.push_back(SPIN::Paired);
        Phi_a.push_back(SPIN::Paired);
        mpi::distribute(Phi_a);

        if (mpi::my_orb(Phi_a[0])) Phi_a[0].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi_a[1])) Phi_a[1].alloc(NUMBER::Imag);
        if (mpi::my_orb(Phi_a[0])) mrcpp::project<3>(prec, Phi_a[0].real(), f1);
        if (mpi::my_orb(Phi_a[1])) mrcpp::project<3>(prec, Phi_a[1].imag(), f2);

        OrbitalVector Phi_b;
        Phi_b.push_back(SPIN::Paired);
        Phi_b.push_back(SPIN::Paired);
        mpi::distribute(Phi_b);

        if (mpi::my_orb(Phi_b[0])) Phi_b[0].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi_b[1])) Phi_b[1].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi_b[0])) mrcpp::project<3>(prec, Phi_b[0].real(), f3);
        if (mpi::my_orb(Phi_b[1])) mrcpp::project<3>(prec, Phi_b[1].real(), f4);

        OrbitalVector Phi_c = add(1.0, Phi_a, -1.0, Phi_b);

        ComplexVector int_a = get_integrals(Phi_a);
        ComplexVector int_b = get_integrals(Phi_b);
        ComplexVector int_c = get_integrals(Phi_c);

        REQUIRE(int_c(0).real() == Approx(int_a[0].real() - int_b[0].real()));
        REQUIRE(int_c(0).imag() == Approx(int_a[0].imag() - int_b[0].imag()));
        REQUIRE(int_c(1).real() == Approx(int_a[1].real() - int_b[1].real()));
        REQUIRE(int_c(1).imag() == Approx(int_a[1].imag() - int_b[1].imag()));

        free(Phi_a);
        free(Phi_b);
        free(Phi_c);
    }

    SECTION("orbital transformations") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        mpi::distribute(Phi);

        if (mpi::my_orb(Phi[0])) {
            Phi[0].alloc(NUMBER::Total);
            mrcpp::project<3>(prec, Phi[0].real(), f1);
            mrcpp::project<3>(prec, Phi[0].imag(), f2);
        }

        if (mpi::my_orb(Phi[1])) {
            Phi[1].alloc(NUMBER::Total);
            mrcpp::project<3>(prec, Phi[1].real(), f3);
            mrcpp::project<3>(prec, Phi[1].imag(), f4);
        }

        orthogonalize(Phi);
        normalize(Phi);

        double theta = 0.5;
        ComplexMatrix U(Phi.size(), Phi.size());
        U(0, 0) = cos(theta);
        U(0, 1) = -sin(theta);
        U(1, 0) = sin(theta);
        U(1, 1) = cos(theta);

        SECTION("unitary transformation") {
            OrbitalVector Psi = rotate(U, Phi, prec);
            ComplexMatrix S = orbital::calc_overlap_matrix(Psi);
            for (int i = 0; i < S.rows(); i++) {
                for (int j = 0; j < S.cols(); j++) {
                    if (i == j) REQUIRE(std::abs(S(i, j)) == Approx(1.0));
                    if (i != j) REQUIRE(std::abs(S(i, j)) < thrs);
                }
            }
            free(Psi);
        }

        free(Phi);
    }
}

} //namespace orbital_vector_tests
