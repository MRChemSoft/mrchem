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

#include "tensor/RankOneOperator.h"

#include "qmoperators/QMDerivative.h"
#include "qmoperators/one_electron/NablaOperator.h"

namespace mrchem {

/** @class HessianOperator
 *
 * @brief Hessian operator
 *
 * This Hessian operator stores the second derivatives of Orbitals in an OrbitalVector. The ordering of the indices is
 * [xx, yy, zz, yz, xz, xy].
 */
class HessianOperator final : public RankOneOperator<6> {
public:
    HessianOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D1, std::shared_ptr<mrcpp::DerivativeOperator<3>> D2, double prec) {
        bool imag = false;
        NablaOperator N1 = NablaOperator(D1, imag);
        NablaOperator N2 = NablaOperator(D2, imag);
        N1.setup(prec);
        N2.setup(prec);

        // Invoke operator= to assign *this operator
        RankOneOperator<6> &d = (*this);
        d[0] = N2[0];
        d[1] = N2[1];
        d[2] = N2[2];
        d[3] = N1[1] * N1[2];
        d[4] = N1[0] * N1[2];
        d[5] = N1[0] * N1[1];

        d[0].name() = "del[x]del[x]";
        d[1].name() = "del[y]del[y]";
        d[2].name() = "del[z]del[z]";
        d[3].name() = "del[y]del[z]";
        d[4].name() = "del[x]del[z]";
        d[5].name() = "del[x]del[y]";
    }
};

} // namespace mrchem