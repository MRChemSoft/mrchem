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

#include "H_E_dip.h"

/** @class ElectricFieldOperator
 *
 * @brief External electric field operator
 *
 * An external electric field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the electric field vector.
 *
 * It implements also explicit trace functions for the nuclear contributions
 *
 */

namespace mrchem {

/** @brief constructor
 *
 * @param[in] f the external electric field
 *
 */
class ElectricFieldOperator final : public RankZeroOperator {
public:
    ElectricFieldOperator(const Eigen::Vector3d &f, const mrcpp::Coord<3> &o)
            : ElectricFieldOperator(std::array<double, 3>{f[0], f[1], f[2]}, H_E_dip(o)) {}

    ElectricFieldOperator(const Eigen::Vector3d &f, H_E_dip mu)
            : ElectricFieldOperator(std::array<double, 3>{f[0], f[1], f[2]}, mu) {}

    ElectricFieldOperator(const std::array<double, 3> &f, const mrcpp::Coord<3> &o)
            : ElectricFieldOperator(f, H_E_dip(o)) {}

    ElectricFieldOperator(const std::array<double, 3> &f, H_E_dip mu) {
        // Invoke operator= to assign *this operator
        RankZeroOperator &HEF = (*this);
        HEF = -f[0] * mu[0] - f[1] * mu[1] - f[2] * mu[2];
        HEF.name() = "E . mu_E";
    }
};

} // namespace mrchem
