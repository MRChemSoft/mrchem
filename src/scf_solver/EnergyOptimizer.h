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

#pragma once

#include "GroundStateSolver.h"

/** @class EnergyOptimizer
 *
 * @brief Ground state SCF optimization without kinetic energy operator
 *
 * This ground state SCF solver completely avoids the use of the kinetic energy
 * operator in the solution algorithm, in particular by computing a potential
 * update to the Fock matrix at the given iteration. This algorithm is considerably
 * less robust than the regular OrbitalOptimizer since there is no KAIN accelerator,
 * and is meant to be used only in the final iterations of the SCF procedure in order
 * to ensure quadratic precision of the energy relative to the orbitals (the effect
 * of this is currently questionable, after the improved implementation of derivative
 * operators in MRCPP).
 */

namespace mrchem {

class EnergyOptimizer final : public GroundStateSolver {
public:
    EnergyOptimizer(HelmholtzVector &h);

    void setup(FockOperator &fock,
               OrbitalVector &phi,
               ComplexMatrix &F,
               FockOperator &fock_np1,
               OrbitalVector &phi_np1);
    void clear();

    bool optimize();

protected:
    FockOperator *fOper_np1;     ///< Next iteration Fock operator (pointer to external object)
    OrbitalVector *orbitals_np1; ///< Next iteration orbitals (pointer to external object)

    ComplexMatrix calcFockMatrixUpdate(double prec, OrbitalVector &dPhi_n);
};

} // namespace mrchem
