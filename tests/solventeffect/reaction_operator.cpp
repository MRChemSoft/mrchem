/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <array>
#include <cmath>
#include <vector>

#include "MRCPP/MWFunctions"
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

namespace reaction_operator {

TEST_CASE("ReactionOperator", "[reaction_operator]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-8;

    // initialize operators, kain history and linearity
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    auto history = 4;
    // Initialize cavity and dielectric constants
    std::vector<mrcpp::Coord<3>> coords = {{0.0, 0.0, 0.0}};
    std::vector<double> R = {1.0};
    double slope = 0.2;
    auto sphere = std::make_shared<Cavity>(coords, R, slope);

    auto PT = std::make_shared<PeriodicTable>();
    auto N = std::make_shared<Nucleus>(PT->getElement(1), coords[0]);
    auto molecule = std::make_shared<Nuclei>();
    molecule->push_back(*N);

    auto Phi_p = std::make_shared<OrbitalVector>();

    Phi_p->push_back(Orbital(SPIN::Paired));
    HydrogenFunction f(1, 0, 0);
    qmfunction::project((*Phi_p)[0], f, NUMBER::Real, prec);

    ReactionOperator Reo(P_p, D_p, sphere, *molecule, Phi_p, history, 1.0, 2.0, false);
    Reo.setup(prec);
    auto gamma_test = Reo.getGamma();
    REQUIRE(Reo.getTotalEnergy() == Approx(-0.210043571463).epsilon(thrs));
    REQUIRE(gamma_test.integrate().real() == Approx(-0.18156088037).epsilon(thrs));
    Reo.clear();
}
} // namespace reaction_operator
