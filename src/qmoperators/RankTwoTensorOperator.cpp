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

#include "RankTwoTensorOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

namespace mrchem {

template<int I, int J>
ComplexMatrix RankTwoTensorOperator<I,J>::operator()(Orbital bra, Orbital ket) {
    RankTwoTensorOperator<I,J> &O = *this;
    ComplexMatrix out(I,J);
    for (int i = 0; i < I; i++) {
        out.row(i) = O[i](bra, ket);
    }
    return out;
}

template<int I, int J>
ComplexMatrix RankTwoTensorOperator<I,J>::trace(OrbitalVector &phi) {
    RankTwoTensorOperator<I,J> &O = *this;
    ComplexMatrix out(I,J);
    for (int i = 0; i < I; i++) {
        out.row(i) = O[i].trace(phi);
    }
    return out;
}

template<int I, int J>
ComplexMatrix RankTwoTensorOperator<I,J>::trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
    RankTwoTensorOperator<I,J> &O = *this;
    ComplexMatrix out(I,J);
    for (int i = 0; i < I; i++) {
        out.row(i) = O[i].trace(phi, x, y);
    }
    return out;
}

} //namespace mrchem

template class mrchem::RankTwoTensorOperator<3,3>;
