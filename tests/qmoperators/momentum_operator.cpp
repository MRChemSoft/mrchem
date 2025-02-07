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

#include "analyticfunctions/HarmonicOscillatorFunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/MomentumOperator.h"

using namespace mrchem;
using namespace orbital;

namespace momentum_operator {

TEST_CASE("MomentumOperator", "[momentum_operator]") {
    const double prec = 1.0e-3;
    const double thrs = prec * prec;

    int nFuncs = 3;
    OrbitalVector Phi;
    for (int n = 0; n < nFuncs; n++) Phi.push_back(Orbital(SPIN::Paired));
    Phi.distribute();

    for (int n = 0; n < nFuncs; n++) {
        int nu[3] = {n, 0, 0};
        HarmonicOscillatorFunction f(nu);
        if (mrcpp::mpi::my_orb(Phi[n])) mrcpp::cplxfunc::project(Phi[n], f, NUMBER::Real, prec);
    }

    // reference values for harmonic oscillator eigenfunctions
    DoubleMatrix ref(nFuncs, nFuncs);
    for (int i = 0; i < nFuncs; i++) {
        for (int j = 0; j < nFuncs; j++) {
            ref(i, j) = 0.0;
            if (i == j + 1) ref(i, j) = std::sqrt(i / 2.0);
            if (i == j - 1) ref(i, j) = -std::sqrt(j / 2.0);
        }
    }
    auto D = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
    MomentumOperator p(D);
    p.setup(prec);
    SECTION("vector apply") {
        OrbitalVector xPhi = p[0](Phi);
        ComplexMatrix X = orbital::calc_overlap_matrix(Phi, xPhi);
        for (int i = 0; i < X.rows(); i++) {
            for (int j = 0; j < X.cols(); j++) { REQUIRE(std::abs(X(i, j).imag() - ref(i, j)) < thrs); }
        }
    }
    SECTION("expectation matrix ") {
        ComplexMatrix X = p[0](Phi, Phi);
        for (int i = 0; i < X.rows(); i++) {
            for (int j = 0; j < X.cols(); j++) { REQUIRE(std::abs(X(i, j).imag() - ref(i, j)) < thrs); }
        }
    }
    p.clear();
}

} // namespace momentum_operator
