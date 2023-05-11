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

#include "RankZeroOperator.h"
#include "TensorOperator.h"

namespace mrchem {

using Orbital = mrcpp::ComplexFunction;

/** @class RankOneOperator
 *
 *  @brief Vector of RankZeroOperator
 *
 * This class provides a base for all vector operators, and implements some simple
 * collective operations returning vector quantities.
 *
 */

template <int I> class RankOneOperator : public TensorOperator<I, RankZeroOperator> {
public:
    RankOneOperator<I> operator()(RankZeroOperator B);
    OrbitalVector operator()(Orbital phi);
    ComplexVector operator()(Orbital bra, Orbital ket);
    ComplexVector trace(OrbitalVector &phi);
    ComplexVector trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y);
    ComplexVector trace(const Nuclei &nucs);
};

namespace tensor {
template <int I> RankZeroOperator dot(RankOneOperator<I> A, RankOneOperator<I> B);
RankOneOperator<3> cross(RankOneOperator<3> A, RankOneOperator<3> B);
} // namespace tensor

} // namespace mrchem
