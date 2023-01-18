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

#pragma once

#include "qmoperators/two_electron/XCPotential.h"

/** @class XCPotential
 *
 * @brief Exchange-Correlation potential defined by a particular (spin) density
 *
 * The XC potential is computed by mapping of the density through a XC functional,
 * provided by the XCFun library. There are two ways of defining the density:
 *
 *  1) Use getDensity() prior to setup() and build the density as you like.
 *  2) Provide a default set of orbitals in the constructor that is used to
 *     compute the density on-the-fly in setup().
 *
 * If a set of orbitals has NOT been given in the constructor, the density
 * MUST be explicitly computed prior to setup(). The density will be computed
 * on-the-fly in setup() ONLY if it is not already available. After setup() the
 * operator will be fixed until clear(), which deletes both the density and the
 * potential.
 *
 * LDA and GGA functionals are supported as well as two different ways to compute
 * the XC potentials: either with explicit derivatives or gamma-type derivatives.
 */

namespace mrchem {

class XCPotentialD2 final : public XCPotential {
public:
    XCPotentialD2(std::unique_ptr<mrdft::MRDFT> &F, std::shared_ptr<OrbitalVector> Phi, std::shared_ptr<OrbitalVector> X, std::shared_ptr<OrbitalVector> Y, bool mpi_shared = false);

private:
    std::shared_ptr<OrbitalVector> orbitals_x; ///< 1st external set of perturbed orbitals used to build the density
    std::shared_ptr<OrbitalVector> orbitals_y; ///< 2nd external set of perturbed orbitals used to build the density

    mrcpp::FunctionTreeVector<3> setupDensities(double prec, mrcpp::FunctionTree<3> &grid);
};

} // namespace mrchem
