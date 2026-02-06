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

#include "tensor/RankZeroOperator.h"

#include "H_B_dip.h"

/** @class MagneticFieldOperatorImag
 *
 * @brief Imaginary component of the external magnetic field operator
 *
 * An external magnetic field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the magnetic field vector.
 *
 */

namespace mrchem {

class MagneticFieldOperatorImag final : public RankZeroOperator {
public:
    MagneticFieldOperatorImag(const Eigen::Vector3d &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : MagneticFieldOperatorImag(std::array<double, 3>{f[0], f[1], f[2]}, H_B_dip(D, o)) {}

    MagneticFieldOperatorImag(const Eigen::Vector3d &f, H_B_dip mu)
            : MagneticFieldOperatorImag(std::array<double, 3>{f[0], f[1], f[2]}, mu) {}

    MagneticFieldOperatorImag(const std::array<double, 3> &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : MagneticFieldOperatorImag(f, H_B_dip(D, o)) {}

    MagneticFieldOperatorImag(const std::array<double, 3> &f, H_B_dip mu) {
        // -i * B * L_0
        RankZeroOperator MF_dip = f[0] * mu[0] + f[1] * mu[1] + f[2] * mu[2];
        MF_dip.setImag(true);

        // Invoke operator= to assign *this operator
        RankZeroOperator &HMF = (*this);
        //HMF = f[0] * mu[0] + f[1] * mu[1] + f[2] * mu[2];
        HMF = MF_dip;
        HMF.name() = "B . mu_B";
    }
};

} // namespace mrchem
