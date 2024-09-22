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

#include "RankOneOperator.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

namespace mrchem {

template <int I> RankOneOperator<I> RankOneOperator<I>::operator()(RankZeroOperator B) {
    RankOneOperator<I> &A = *this;
    RankOneOperator<I> out;
    for (int i = 0; i < I; i++) out[i] = A[i](B);
    return out;
}

template <int I> OrbitalVector RankOneOperator<I>::operator()(Orbital phi) {
    RankOneOperator<I> &O = *this;
    OrbitalVector out;
    for (int i = 0; i < I; i++) out.push_back(O[i](phi));
    return out;
}

template <int I> ComplexVector RankOneOperator<I>::operator()(Orbital bra, Orbital ket) {
    RankOneOperator<I> &O = *this;
    ComplexVector out(I);
    for (int i = 0; i < I; i++) out(i) = O[i](bra, ket);
    return out;
}

template <int I> ComplexVector RankOneOperator<I>::trace(OrbitalVector &phi) {
    RankOneOperator<I> &O = *this;
    ComplexVector out = ComplexVector::Zero(I);
    for (int i = 0; i < I; i++) out(i) = O[i].trace(phi);
    return out;
}

template <int I> ComplexVector RankOneOperator<I>::trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
    RankOneOperator<I> &O = *this;
    ComplexVector out = ComplexVector::Zero(I);
    for (int i = 0; i < I; i++) out(i) = O[i].trace(phi, x, y);
    return out;
}

template <int I> ComplexVector RankOneOperator<I>::trace(const Nuclei &nucs) {
    RankOneOperator<I> &O = *this;
    ComplexVector out = ComplexVector::Zero(I);
    for (int i = 0; i < I; i++) out(i) = O[i].trace(nucs);
    return out;
}

} // namespace mrchem

template class mrchem::RankOneOperator<3>;
template class mrchem::RankOneOperator<6>;
