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

#include "catch.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "MRCPP/MWFunctions"
#include <MRCPP/MWOperators>

#include "mrchem.h"
#include "parallel.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "chemistry/Element.h"
#include "chemistry/Nucleus.h"
#include "chemistry/PeriodicTable.h"
#include "environment/Cavity.h"
#include "environment/Permittivity.h"
#include "environment/SCRF.h"
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

    // initialize operators
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    std::shared_ptr<mrcpp::DerivativeOperator<3>> D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);

    // initialize spherical cavity
    double slope = 0.2;
    std::vector<double> radius = {1.0};
    std::vector<mrcpp::Coord<3>> coords = {{0.0, 0.0, 0.0}};
    Cavity sphere(coords, radius, slope);

    // initialize dielectric function
    double eps_in = 1.0;
    double eps_out = 2.0;
    Permittivity dielectric_func(sphere, eps_in, eps_out, "exponential");

    // initialize molecule containing single hydrogen
    PeriodicTable PT;
    Nucleus H(PT.getElement(1), coords[0]);
    Nuclei molecule;
    molecule.push_back(H);

    // initialize orbital vector
    auto Phi_p = std::make_shared<OrbitalVector>();
    auto &Phi = *Phi_p;
    Phi.push_back(Orbital(SPIN::Paired));
    mpi::distribute(Phi);

    // project analytic 1s orbital
    HydrogenFunction f(1, 0, 0);
    if (mpi::my_orb(Phi[0])) qmfunction::project(Phi[0], f, NUMBER::Real, prec);

    int kain = 4;
    auto scrf_p = std::make_unique<SCRF>(dielectric_func, molecule, P_p, D_p, prec, kain, 100, true, false, "total");

    auto Reo = std::make_shared<ReactionOperator>(std::move(scrf_p), Phi_p);
    Reo->setup(prec);
    double total_energy = Reo->getTotalEnergy();
    Reo->clear();
    REQUIRE(total_energy == Approx(-0.191434124263).epsilon(thrs));
}

} // namespace reaction_operator
