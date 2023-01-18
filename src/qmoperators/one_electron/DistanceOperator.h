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

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/QMPotential.h"

namespace mrchem {

class DistanceOperator final : public RankZeroOperator {
public:
    /*! @brief DistanceOperator represents the function: 1.0/|r - o|^p
     *  @param p: Power in denominator
     *  @param o: Coordinate of origin
     *  @param proj_prec: Precision for projection of analytic function
     *  @param smooth_prec: Precision for smoothing of analytic function
     */
    DistanceOperator(double p, const mrcpp::Coord<3> &o, double proj_prec, double smooth_prec = -1.0) {
        if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
        if (smooth_prec < 0.0) smooth_prec = proj_prec;

        // Define analytic smoothed 1/r
        double c = detail::nuclear_potential_smoothing(smooth_prec, 1.0);
        NuclearFunction nuc_func;
        nuc_func.push_back("H", o, std::pow(c, p));

        // Raise analytic 1/r to power
        auto f = [p, nuc_func](const mrcpp::Coord<3> &r) -> double {
            double f_r = nuc_func.evalf(r);
            return std::pow(f_r, p);
        };

        // Project analytic potential, building grid for 1/r
        auto r_pow = std::make_shared<QMPotential>(1);
        r_pow->alloc(NUMBER::Real);
        mrcpp::build_grid(r_pow->real(), nuc_func);
        mrcpp::project<3>(proj_prec, r_pow->real(), f);

        std::stringstream o_name;
        o_name << "r^{" << std::setprecision(1) << std::fixed << p << "}";

        // Invoke operator= to assign *this operator
        RankZeroOperator &O = (*this);
        O = r_pow;
        O.name() = o_name.str();
    }
};

} // namespace mrchem
