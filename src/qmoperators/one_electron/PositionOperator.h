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

#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"

namespace mrchem {

class PositionOperator : public RankOneOperator<3> {
public:
    /*! @brief PositionOperator represents the vector operator: hat{r} = [r_x - o_x, r_y - o_y, r_z - o_z]
     *  @param o: Coordinate of origin
     */
    PositionOperator(const mrcpp::Coord<3> &o = {0.0, 0.0, 0.0}) {
        // Define analytic potential
        auto f_x = [o](const mrcpp::Coord<3> &r) -> double { return (r[0] - o[0]); };
        auto f_y = [o](const mrcpp::Coord<3> &r) -> double { return (r[1] - o[1]); };
        auto f_z = [o](const mrcpp::Coord<3> &r) -> double { return (r[2] - o[2]); };

        auto r_x = std::make_shared<QMPotential>(1);
        auto r_y = std::make_shared<QMPotential>(1);
        auto r_z = std::make_shared<QMPotential>(1);

        // Project analytic potential (exact on root scale, thus no prec)
        qmfunction::project(*r_x, f_x, NUMBER::Real, -1.0);
        qmfunction::project(*r_y, f_y, NUMBER::Real, -1.0);
        qmfunction::project(*r_z, f_z, NUMBER::Real, -1.0);

        // Invoke operator= to assign *this operator
        RankOneOperator &r = (*this);
        r[0] = r_x;
        r[1] = r_y;
        r[2] = r_z;
        r[0].name() = "r[x]";
        r[1].name() = "r[y]";
        r[2].name() = "r[z]";
    }
};

} // namespace mrchem
