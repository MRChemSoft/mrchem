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

#include "qmoperators/QMSpin.h"

namespace mrchem {

class SpinOperator final : public RankOneOperator<3> {
public:
    SpinOperator() {
        auto s_x = std::make_shared<QMSpin>(0);
        auto s_y = std::make_shared<QMSpin>(1);
        auto s_z = std::make_shared<QMSpin>(2);

        // Invoke operator= to assign *this operator
        RankOneOperator &s = *this;
        s[0] = s_x;
        s[1] = s_y;
        s[2] = s_z;
        s[0].name() = "s[x]";
        s[1].name() = "s[y]";
        s[2].name() = "s[z]";
    }
};

} // namespace mrchem
