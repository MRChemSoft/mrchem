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

namespace mrchem {

/** @class TensorOperator
 *
 *  @brief Placeholder class for QMOperators
 *
 * This class provides a tensor notation for general QM operators. It is the base
 * for RankOneTensorOperator and RankTwoTensorOperator (can easily be extended to
 * higher ranks), and allows the components of the tensor operator to be accessed
 * in the usual way (x,y,z = integers):
 *
 * Rank 1: h[x]
 * Rank 2: h[x][y]
 * Rank 3: h[x][y][z]
 *
 * This is achieved by recursion: RankOneTensorOperator takes RankZeroTensorOperator
 * as second template, RankTwoTensorOperator takes RankOneTensorOperator as second
 * template, etc. (first argument is the dimension of the current rank). Note that
 * RankZeroTensorOperator is NOT derived from this class, but is a more general
 * operator expansion (see also @file qmoperator_fwd.h).
 *
 */

template<int I, class T>
class TensorOperator {
public:
    TensorOperator() { }
    virtual ~TensorOperator() { }

    void setup(double prec) { for (int i = 0; i < I; i++) this->oper[i].setup(prec); }
    void clear() { for (int i = 0; i < I; i++) this->oper[i].clear(); }

    T& operator[](int i) { return this->oper[i]; }
    const T& operator[](int i) const { return this->oper[i]; }

protected:
    T oper[I];
};

} //namespace mrchem
