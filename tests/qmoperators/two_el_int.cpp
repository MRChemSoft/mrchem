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

#include "MRCPP/MWOperators"

#include "mrchem.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/two_electron/two_electron_utils.h"

using namespace mrchem;
using namespace orbital;

TEST_CASE("TwoElInt", "[two_el_int]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-8;

    const int nShells = 2;
    std::vector<int> ns;
    std::vector<int> ls;
    std::vector<int> ms;

    auto Phi_p = std::make_shared<OrbitalVector>();

    OrbitalVector &Phi = *Phi_p;
    ns.push_back(1);
    ls.push_back(0);
    ms.push_back(0);
    Phi.push_back(Orbital(SPIN::Paired));

    ns.push_back(2);
    ls.push_back(0);
    ms.push_back(0);
    Phi.push_back(Orbital(SPIN::Paired));
    Phi.distribute();

    for (int i = 0; i < Phi.size(); i++) {
        HydrogenFunction f(ns[i], ls[i], ms[i]);
        if (mrcpp::mpi::my_func(Phi[i])) mrcpp::project(Phi[i], f, prec);
    }

    Eigen::Tensor<std::complex<double>, 4>two_el_int(nShells, nShells, nShells, nShells);

    double t0000 = 0.6249989698;
    double t1000 = 0.0893550109;
    double t1100 = 0.0219478674;
    double t1010 = 0.2099575786;
    double t1110 = 0.0085813860;
    double t1111 = 0.1503904996;

    SECTION("calc_2elintegrals") {
        two_el_int = calc_2elintegrals(prec, Phi);
        REQUIRE(two_el_int(0,0,0,0).real() == Catch::Approx(t0000).margin(prec));
        REQUIRE(two_el_int(0,0,0,0).imag() < thrs);
        REQUIRE(two_el_int(1,0,0,0).real() == Catch::Approx(t1000).margin(prec));
        REQUIRE(two_el_int(0,1,0,0).real() == Catch::Approx(t1000).margin(prec));
        REQUIRE(two_el_int(0,0,1,0).real() == Catch::Approx(t1000).margin(prec));
        REQUIRE(two_el_int(0,0,0,1).real() == Catch::Approx(t1000).margin(prec));
        REQUIRE(two_el_int(1,1,0,0).real() == Catch::Approx(t1100).margin(prec));
        REQUIRE(two_el_int(1,0,1,0).real() == Catch::Approx(t1010).margin(prec));
        REQUIRE(two_el_int(1,1,1,0).real() == Catch::Approx(t1110).margin(prec));
        REQUIRE(two_el_int(1,1,1,1).real() == Catch::Approx(t1111).margin(prec));
   }

}
