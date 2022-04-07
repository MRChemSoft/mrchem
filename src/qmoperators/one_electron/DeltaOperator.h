/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <MRCPP/Gaussians>

#include "tensor/RankZeroOperator.h"

#include "chemistry/PhysicalConstants.h"

#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"

namespace mrchem {

class DeltaOperator final : public RankZeroOperator {
public:
    /*! @brief DeltaOperator represents the Dirac delta function: delta|r - o|
     *  @param o: Coordinate of origin
     *  @param proj_prec: Precision for projection of analytic function
     *  @param smooth_prec: Precision for smoothing of analytic function
     */
    DeltaOperator(const mrcpp::Coord<3> &o, double proj_prec, double smooth_prec = -1.0) {
        if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
        if (smooth_prec < 0.0) smooth_prec = proj_prec;

        // Define analytic potential
        double beta = 1.0 / smooth_prec;
        double alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        mrcpp::GaussFunc<3> f(beta, alpha, o);

        // Project analytic potential
        auto delta = std::make_shared<QMPotential>(1);
        qmfunction::project(*delta, f, NUMBER::Real, proj_prec);

        // Invoke operator= to assign *this operator
        RankZeroOperator &h = (*this);
        h = delta;
        h.name() = "delta";
    }
};

} // namespace mrchem
