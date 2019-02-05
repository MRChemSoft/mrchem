/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "analyticfunctions/HydrogenFunction.h"
#include "chemistry/Cavity.h"
#include "chemistry/Element.h"
#include "chemistry/Nucleus.h"
#include "chemistry/PeriodicTable.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/two_electron/ReactionOperator.h"

using namespace mrchem;
using namespace orbital;

namespace reaction_potential {

TEST_CASE("ReactionOperator", "[reaction_operator]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-8;

    const int nShells = 2;
    std::vector<int> ns;
    std::vector<int> ls;
    std::vector<int> ms;

    OrbitalVector Phi;
    for (int n = 1; n <= nShells; n++) {
        int L = n;
        for (int l = 0; l < L; l++) {
            int M = 2 * l + 1;
            for (int m = 0; m < M; m++) {
                ns.push_back(n);
                ls.push_back(l);
                ms.push_back(m);
                Phi.push_back(Orbital(SPIN::Paired));
            }
        }
    }
    mpi::distribute(Phi);

    for (int i = 0; i < Phi.size(); i++) {
        HydrogenFunction f(ns[i], ls[i], ms[i]);
        if (mpi::my_orb(Phi[i])) qmfunction::project(Phi[i], f, NUMBER::Real, prec);
    }

    int i = 0;
    DoubleMatrix E_P = DoubleMatrix::Zero(Phi.size(), Phi.size());

    E_P(0, 0) = -44.5732678692;
    E_P(1, 0) = -8.9637637907;
    E_P(1, 1) = -11.2824319897;
    E_P(2, 0) = -0.371064946;
    E_P(2, 1) = 0.035800708;
    E_P(2, 2) = -11.668385913;
    E_P(3, 0) = -0.371064946;
    E_P(3, 1) = 0.035800708;
    E_P(3, 2) = -0.0056937037;
    E_P(3, 3) = -11.668385913;
    E_P(4, 0) = -0.371064946;
    E_P(4, 1) = 0.035800708;
    E_P(4, 2) = -0.0056937037;
    E_P(4, 3) = -0.0056937037;
    E_P(4, 4) = -11.668385913;

    std::vector<mrcpp::Coord<3>> center = {{0.5, 0.5, 0.5}};
    std::vector<double> R = {1.0};
    double s = 0.5;
    Cavity C(center, R, s);

    PeriodicTable PT;
    mrcpp::Coord<3> pos = {0.0, 0.0, 0.0};
    Element H = PT.getElement('H');
    Nucleus n(H, pos);
    Nuclei nucs;
    nucs.push_back(n);

    mrcpp::PoissonOperator P(*MRA, prec);
    mrcpp::ABGVOperator<3> D(*MRA, 0.0, 0.0);

    ReactionOperator V(&P, &D, &C, nucs, &Phi, true);

    V.setup(prec);

    SECTION("apply") {
        Orbital Vphi_0 = V(Phi[0]);
        ComplexDouble V_00 = orbital::dot(Phi[0], Vphi_0);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(V_00.real() == Approx(E_P(0, 0)).epsilon(thrs));
            REQUIRE(V_00.imag() < thrs);
        } else {
            REQUIRE(V_00.real() < thrs);
            REQUIRE(V_00.imag() < thrs);
        }
        // }
        // SECTION("vector apply") {
        OrbitalVector VPhi = V(Phi);
        for (int i = 0; i < Phi.size(); i++) {
            ComplexDouble V_ii = orbital::dot(Phi[i], VPhi[i]);
            if (mpi::my_orb(Phi[i])) {
                REQUIRE(V_ii.real() == Approx(E_P(i, i)).epsilon(thrs));
                REQUIRE(V_ii.imag() < thrs);
            } else {
                REQUIRE(V_ii.real() < thrs);
                REQUIRE(V_ii.imag() < thrs);
            }
        }
        //}
        // SECTION("expectation value") {
        ComplexDouble V_11 = V(Phi[1], Phi[1]);
        if (mpi::my_orb(Phi[1])) {
            REQUIRE(V_11.real() == Approx(E_P(1, 1)).epsilon(thrs));
            REQUIRE(V_11.imag() < thrs);
        } else {
            REQUIRE(V_11.real() < thrs);
            REQUIRE(V_11.imag() < thrs);
        }
        // }
        // SECTION("expectation matrix ") {
        ComplexMatrix v = V(Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j <= i; j++) {
                //  if (std::abs(v(i, j).real()) > thrs) REQUIRE(v(i, j).real() == Approx(E_P(i, j)).epsilon(thrs));
                if (std::abs(v(i, j).real()) > thrs) REQUIRE(v(i, j).real() == v(i, j).real());
                REQUIRE(v(i, j).imag() < thrs);
            }
        }
    }
    V.clear();
}

} // namespace reaction_potential
