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

#include "MRCPP/MWOperators"

#include "mrchem.h"
#include "parallel.h"

#include "analyticfunctions/HarmonicOscillatorFunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/one_electron/KineticOperator.h"

using namespace mrchem;
using namespace orbital;

namespace kinetic_operator {

TEST_CASE("KineticOperator", "[kinetic_operator]") {
    const double prec = 1.0e-3;
    const double thrs = prec * prec;

    int nFuncs = 3;
    OrbitalVector Phi;
    for (int n = 0; n < nFuncs; n++) Phi.push_back(Orbital(SPIN::Paired));
    mpi::distribute(Phi);

    for (int n = 0; n < nFuncs; n++) {
        int nu[3] = {n, 0, 0};
        HarmonicOscillatorFunction f(nu);
        if (mpi::my_orb(Phi[n])) qmfunction::project(Phi[n], f, NUMBER::Real, prec);
    }

    // reference values for harmonic oscillator eigenfunctions
    DoubleVector E_K(nFuncs);
    for (int i = 0; i < nFuncs; i++) {
        //energy = (nu + 1/2)
        double E_x = (i + 0.5);
        double E_y = (0 + 0.5);
        double E_z = (0 + 0.5);

        //virial theorem: <E_K> = <E_P> = E/2
        E_K(i) = 0.5 * (E_x + E_y + E_z);
    }

    mrcpp::ABGVOperator<3> D(*MRA, 0.5, 0.5);
    KineticOperator T(D);
    T.setup(prec);
    SECTION("apply") {
        Orbital Tphi_0 = T(Phi[0]);
        ComplexDouble T_00 = orbital::dot(Phi[0], Tphi_0);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(T_00.real() == Approx(E_K(0)));
            REQUIRE(T_00.imag() < thrs);
        } else {
            REQUIRE(T_00.real() < thrs);
            REQUIRE(T_00.imag() < thrs);
        }
    }
    SECTION("vector apply") {
        OrbitalVector TPhi = T(Phi);
        ComplexMatrix t = orbital::calc_overlap_matrix(Phi, TPhi);
        for (int i = 0; i < Phi.size(); i++) {
            ComplexDouble T_ii = orbital::dot(Phi[i], TPhi[i]);
            if (mpi::my_orb(Phi[i])) {
                REQUIRE(T_ii.real() == Approx(E_K(i)));
                REQUIRE(T_ii.imag() < thrs);
            } else {
                REQUIRE(T_ii.real() < thrs);
                REQUIRE(T_ii.imag() < thrs);
            }
        }
    }
    SECTION("expectation value") {
        ComplexDouble T_00 = T(Phi[0], Phi[0]);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(T_00.real() == Approx(E_K(0)));
            REQUIRE(T_00.imag() < thrs);
        } else {
            REQUIRE(T_00.real() < thrs);
            REQUIRE(T_00.imag() < thrs);
        }
    }
    SECTION("expectation matrix ") {
        ComplexMatrix t = T(Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            REQUIRE(t(i, i).real() == Approx(E_K(i)));
            REQUIRE(t(i, i).imag() < thrs);
        }
    }
    T.clear();
}

} //namespace kinetic_operator
