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

#include "TensorOperator.h"
#include "RankOneTensorOperator.h"

namespace mrchem {

class Orbital;

/** @class RankTwoTensorOperator
 *
 *  @brief Matrix of RankZeroTensorOperator
 *
 * This class provides a base for all matrix operators (implemented as a vector
 * of vectors), and implements some simple collective operations returning matrix
 * quantities.
 *
 */

template<int I, int J>
class RankTwoTensorOperator : public TensorOperator<I, RankOneTensorOperator<J> > {
public:
    ComplexMatrix operator()(Orbital bra, Orbital ket);
    ComplexMatrix trace(OrbitalVector &phi);
    ComplexMatrix trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y);
};

} //namespace mrchem
