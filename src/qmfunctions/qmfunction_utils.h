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

#include "mrchem.h"
#include "qmfunction_fwd.h"

namespace mrchem {
namespace qmfunction {

ComplexDouble dot(mrcpp::ComplexFunction bra, mrcpp::ComplexFunction ket);
ComplexDouble node_norm_dot(mrcpp::ComplexFunction bra, mrcpp::ComplexFunction ket, bool exact);
void deep_copy(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction &inp);
void add(mrcpp::ComplexFunction &out, ComplexDouble a, mrcpp::ComplexFunction inp_a, ComplexDouble b, mrcpp::ComplexFunction inp_b, double prec);
void project(mrcpp::ComplexFunction &out, std::function<double(const mrcpp::Coord<3> &r)> f, int type, double prec);
void project(mrcpp::ComplexFunction &out, mrcpp::RepresentableFunction<3> &f, int type, double prec);
void project(mrcpp::ComplexFunction &out, mrcpp::GaussExp<3> &f, int type, double prec);
void multiply(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction inp_a, mrcpp::ComplexFunction inp_b, double prec, bool absPrec = false, bool useMaxNorms = false);
void multiply_real(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction inp_a, mrcpp::ComplexFunction inp_b, double prec, bool absPrec = false, bool useMaxNorms = false);
void multiply_imag(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction inp_a, mrcpp::ComplexFunction inp_b, double prec, bool absPrec = false, bool useMaxNorms = false);
void linear_combination(mrcpp::ComplexFunction &out, const ComplexVector &c, ComplexFunctionVector &inp, double prec);

} // namespace qmfunction
} // namespace mrchem
