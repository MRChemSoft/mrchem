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

#include "tensor_utils.h"

#include "RankOneOperator.h"
#include "RankTwoOperator.h"
#include "RankZeroOperator.h"

namespace mrchem {

RankOneOperator<3> tensor::cross(RankOneOperator<3> A, RankOneOperator<3> B) {
    RankOneOperator<3> out;
    out[0] = A[1](B[2]) - A[2](B[1]);
    out[1] = A[2](B[0]) - A[0](B[2]);
    out[2] = A[0](B[1]) - A[1](B[0]);
    return out;
}

template <int I> RankZeroOperator tensor::dot(RankOneOperator<I> A, RankOneOperator<I> B) {
    RankZeroOperator out;
    for (int i = 0; i < I; i++) out += A[i](B[i]);
    return out;
}

template <int I, int J> RankTwoOperator<I, J> tensor::outer(RankOneOperator<I> A, RankOneOperator<J> B) {
    RankTwoOperator<I, J> out;
    for (int i = 0; i < I; i++)
        for (int j = 0; j < J; j++) out[i][j] = A[i](B[j]);
    return out;
}

// /*
//  * @brief: shuffles the indices of a spinor, simulating the application of a Dirac matrix to it
//  * alpha represents the index of the Dirac matrices.
//  * For scalar operators, it is unused.
//  * For 2 component (Weyl/Pauli) spinors, alpha = 0,1,2,3 corresponds to indentiy, sigma_x, y and z respectively.
//  * For 4 component (Dirac) spinors, alpha = 0,1,2,3 corresponds to indentiy, alpha_x, y and z respectively,
//  * 
// */
// int shuffle_Dirac(int alpha, int j) {
//     if (alpha > 0) { //not identity
        
//     }
//     else if (alpha == 3) { //gamma 3 (alpha_z) or sigma_z depending on the number of spinor components
 
//     }
//     return j;
// }

namespace tensor {
template RankZeroOperator dot<3>(RankOneOperator<3> A, RankOneOperator<3> B);
template RankTwoOperator<3, 3> outer<3, 3>(RankOneOperator<3> A, RankOneOperator<3> B);
} // namespace tensor

} // namespace mrchem
