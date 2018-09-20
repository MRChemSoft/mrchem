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

#include "SCF.h"

/** @class LinearResponseSolver
 *
 */

namespace mrchem {

class LinearResponseSolver final : public SCF {
public:
    LinearResponseSolver(HelmholtzVector &h, Accelerator *k_x = nullptr, Accelerator *k_y = nullptr);

    void setupUnperturbed(double prec, FockOperator *fock, OrbitalVector *Phi, ComplexMatrix *F);
    void clearUnperturbed();

    void setup(FockOperator *fock, OrbitalVector *X);
    void setup(double omega, FockOperator *fock, OrbitalVector *X, OrbitalVector *Y);
    void clear();

    bool optimize();

protected:
    bool dynamic;
    double frequency;

    FockOperator *fOper_0;
    FockOperator *fOper_1;

    ComplexMatrix *fMat_0;
    ComplexMatrix *fMat_x;
    ComplexMatrix *fMat_y;

    OrbitalVector *orbitals_0;
    OrbitalVector *orbitals_x;
    OrbitalVector *orbitals_y;

    Accelerator *kain_x;
    Accelerator *kain_y;

    // clang-format off
    OrbitalVector setupHelmholtzArguments(OrbitalVector &dPhi,
                                          const ComplexMatrix &M,
                                          bool adjoint);
    // clang-format on

    void printProperty() const;
    double calcProperty();
};

} // namespace mrchem
