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

#include "analyticfunctions/HydrogenFunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"

using namespace mrchem;
using namespace orbital;
using namespace std;

namespace coulomb_hessian {

TEST_CASE("CoulombHessian", "[coulomb_hessian]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-8;

    std::vector<int> ns;
    std::vector<int> ls;
    std::vector<int> ms;
    
    OrbitalVector Phi;
    ns.push_back(1);
    ls.push_back(0);
    ms.push_back(0);
    Phi.push_back(Orbital(SPIN::Paired));

    ns.push_back(2);
    ls.push_back(0);
    ms.push_back(0);
    Phi.push_back(Orbital(SPIN::Paired));

    mpi::distribute(Phi);

    for (int i = 0; i < Phi.size(); i++) {
        HydrogenFunction f(ns[i], ls[i], ms[i]);
        if (mpi::my_orb(Phi[i])) qmfunction::project(Phi[i], f, NUMBER::Real, prec);
    }

    std::vector<int> ns_x;
    std::vector<int> ls_x;
    std::vector<int> ms_x;

    OrbitalVector Phi_x;
    ns_x.push_back(2);
    ls_x.push_back(0);
    ms_x.push_back(0);
    Phi_x.push_back(Orbital(SPIN::Paired));

    ns_x.push_back(2);
    ls_x.push_back(1);
    ms_x.push_back(1);
    Phi_x.push_back(Orbital(SPIN::Paired));

    mpi::distribute(Phi_x);

    for (int i = 0; i < Phi_x.size(); i++) {
        HydrogenFunction f(ns_x[i], ls_x[i], ms_x[i]);
        if (mpi::my_orb(Phi_x[i])) qmfunction::project(Phi_x[i], f, NUMBER::Real, prec);
    }

    int i = 0;
    DoubleMatrix E_P = DoubleMatrix::Zero(Phi.size(), Phi.size());

    E_P(0,0) = 0.356407048592;
    E_P(0,1) = 0.0873913761675;
    E_P(1,0) = 0.0873913761675;
    E_P(1,1) = 0.0341665770117;

    mrcpp::PoissonOperator P(*MRA, prec);
    CoulombOperator V(&P, &Phi, &Phi_x, &Phi_x);

    V.setup(prec);
    SECTION("apply") {
        Orbital Vphi_0 = V(Phi[0]);
        ComplexDouble V_00 = orbital::dot(Phi[0], Vphi_0);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(V_00.real() == Approx(E_P(0,0)).epsilon(thrs));
            REQUIRE(V_00.imag() < thrs);
        } else {
            REQUIRE(V_00.real() < thrs);
            REQUIRE(V_00.imag() < thrs);
        }
    }
    SECTION("vector apply") {
        OrbitalVector VPhi = V(Phi);
        for (int i = 0; i < Phi.size(); i++) {
            ComplexDouble V_ii = orbital::dot(Phi[i], VPhi[i]);
            if (mpi::my_orb(Phi[i])) {
                REQUIRE(V_ii.real() == Approx(E_P(i,i)).epsilon(thrs));
                REQUIRE(V_ii.imag() < thrs);
            } else {
                REQUIRE(V_ii.real() < thrs);
                REQUIRE(V_ii.imag() < thrs);
            }
        }
    }
    SECTION("expectation value") {
        ComplexDouble V_00 = V(Phi[0], Phi[0]);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(V_00.real() == Approx(E_P(0,0)).epsilon(thrs));
            REQUIRE(V_00.imag() < thrs);
        } else {
            REQUIRE(V_00.real() < thrs);
            REQUIRE(V_00.imag() < thrs);
        }
    }
    SECTION("expectation matrix ") {
        ComplexMatrix v = V(Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j <= i; j++) {
                if (std::abs(v(i, j).real()) > thrs) REQUIRE(v(i, j).real() == Approx(E_P(i,j)).epsilon(thrs));
                REQUIRE(v(i, j).imag() < thrs);
            }
        }
    }
    V.clear();
}

} // namespace nuclear_potential
