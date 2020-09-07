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

#pragma once

#include "Accelerator.h"

/** @class KAIN
 *
 * This class implements the Krylov subspace Accelerated Inexact Newton (KAIN)
 * method as described by R.J. Harrison (J. Comput. Chem. 25, 328, 2004).
 *
 * The SCF problem is formulated as finding the root of the function
 *
 * \f$ f(x^n) = -2H[x^n] - x^n \f$
 *
 * where \f$ x^n \f$ is the vector of orbitals, possibly appended by the
 * Fock matrix \f$ x^n = (\phi^n_0, \phi^n_1, \dots \phi^n_N, F^n) \f$
 */

namespace mrchem {

class KAIN final : public Accelerator {
public:
    KAIN(int max, int min = 0, bool sep = false, std::vector<double> scaling = {})
            : Accelerator(max, min, sep)
            , scaling(scaling) {}

protected:
    void setupLinearSystem() override;
    // clang-format off
    void expandSolution(double prec,
                        OrbitalVector &Phi,
                        OrbitalVector &dPhi,
                        ComplexMatrix *F,
                        ComplexMatrix *dF) override;
    // clang-format on
private:
    std::vector<double> scaling;
};

} // namespace mrchem
