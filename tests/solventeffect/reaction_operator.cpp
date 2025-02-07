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

#include <array>
#include <cmath>
#include <vector>

#include "MRCPP/MWFunctions"
#include <MRCPP/MWOperators>

#include "mrchem.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "chemistry/Element.h"
#include "chemistry/Nucleus.h"
#include "chemistry/PeriodicTable.h"
#include "chemistry/chemistry_utils.h"
#include "environment/Cavity.h"
#include "environment/GPESolver.h"
#include "environment/Permittivity.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
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
    auto sphere = std::make_shared<Cavity>(coords, radius, slope);

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
    Phi.distribute();

    // project analytic 1s orbital
    HydrogenFunction f(1, 0, 0);
    if (mrcpp::mpi::my_orb(Phi[0])) mrcpp::cplxfunc::project(Phi[0], f, NUMBER::Real, prec);

    auto rho_nuc = chemistry::compute_nuclear_density(prec, molecule, 100);

    int kain = 4;
    auto scrf_p = std::make_unique<GPESolver>(dielectric_func, rho_nuc, P_p, D_p, kain, 100, false, SCRFDensityType::TOTAL);

    auto Reo = std::make_shared<ReactionOperator>(std::move(scrf_p), Phi_p);
    Reo->setup(prec);

    Density rho_el(false);
    density::compute(prec, rho_el, Phi, DensityType::Total);
    rho_el.rescale(-1.0);

    auto [Er_nuc, Er_el] = Reo->getSolver()->computeEnergies(rho_el);
    auto total_energy = Er_nuc + Er_el;
    Reo->clear();
    REQUIRE(total_energy == Catch::Approx(-1.022729683846e-01).epsilon(thrs));
}

} // namespace reaction_operator
