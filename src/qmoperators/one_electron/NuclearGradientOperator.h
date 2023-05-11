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

#include "analyticfunctions/NuclearFunction.h"
#include "analyticfunctions/NuclearGradientFunction.h"
#include "qmoperators/QMPotential.h"

namespace mrchem {

class NuclearGradientOperator final : public RankOneOperator<3> {
public:
    /*!
     *  @brief NuclearGradientOperator represents the vector potential: Z * hat{r}/|r - o|^3
     *  @param z: Nuclear charge of nucleus
     *  @param o: Coordinate of origin
     *  @param proj_prec: Precision for projection of analytic function
     *  @param c: Smoothing parameter for analytic function
     */
    NuclearGradientOperator(double z, const mrcpp::Coord<3> &o, double proj_prec, double c) {
        if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");

        // Define analytic potential
        NuclearGradientFunction f_x(0, z, o, c);
        NuclearGradientFunction f_y(1, z, o, c);
        NuclearGradientFunction f_z(2, z, o, c);

        auto x_rm3 = std::make_shared<QMPotential>(1);
        auto y_rm3 = std::make_shared<QMPotential>(1);
        auto z_rm3 = std::make_shared<QMPotential>(1);

        // Project analytic potential
        mrcpp::cplxfunc::project(*x_rm3, f_x, NUMBER::Real, proj_prec);
        mrcpp::cplxfunc::project(*y_rm3, f_y, NUMBER::Real, proj_prec);
        mrcpp::cplxfunc::project(*z_rm3, f_z, NUMBER::Real, proj_prec);

        // Invoke operator= to assign *this operator
        RankOneOperator &v = (*this);
        v[0] = x_rm3;
        v[1] = y_rm3;
        v[2] = z_rm3;
        v[0].name() = "x/r^3";
        v[1].name() = "y/r^3";
        v[2].name() = "z/r^3";
    }
};

} // namespace mrchem
