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

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "QMPotential.h"
#include "qmfunctions/Orbital.h"

#include "QMDerivative.h"
#include "QMIdentity.h"
#include "QMSpin.h"
#include "utils/print_utils.h"

using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

using QMOperator_p = std::shared_ptr<mrchem::QMOperator>;

namespace mrchem {

/** @brief constructor
 *
 * @param adap: extra refinement in output
 *
 * Initializes the mrcpp::CompFunction<3> with NULL pointers for both real and imaginary part.
 * These must be computed in setup() of derived classes. The initial output grid
 * in application will be a copy of the input orbital but NOT a copy of the
 * potential grid. The argument sets how many extra refinement levels is allowed
 * beyond this initial refinement.
 */
QMPotential::QMPotential(int adap, bool shared)
        : mrcpp::CompFunction<3>(0, shared)
        , QMOperator()
        , adap_build(adap) {}

QMPotential::QMPotential(const QMPotential &inp)
        : mrcpp::CompFunction<3>(0, inp.isShared())
        , QMOperator()
        , adap_build(inp.adap_build) {}

/** @brief apply potential
 *
 * @param inp: orbital on which to apply
 *
 * Computes a new orbital that is the product of this potential and the input
 * orbital. Orbital parameters are copied from input.
 */
Orbital QMPotential::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out;
    calc(out, inp, false);

    return out;
}

/** @brief apply complex conjugate potential
 *
 * @param inp: orbital on which to apply
 *
 * Computes a new orbital that is the product of the complex conjugate of this
 * potential and the input orbital. Orbital parameters are copied from input.
 */
Orbital QMPotential::dagger(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out = inp.paramCopy();
    calc(out, inp, true);
    return out;
}

QMOperatorVector QMPotential::apply(QMOperator_p &O) {
    QMIdentity *I = dynamic_cast<QMIdentity *>(&(*O));
    QMPotential *V_inp = dynamic_cast<QMPotential *>(&(*O));

    QMOperatorVector out;
    if (I) {
        // O == identity: skip it
        auto V_out = std::make_shared<QMPotential>(*this);
        mrcpp::deep_copy(*V_out, *this);
        out.push_back(V_out);
    } else if (V_inp) {
        // O == potential: merge into single potential
        auto V_out = std::make_shared<QMPotential>(*this);
        calc(*V_out, *V_inp, false);
        out.push_back(V_out);
    } else {
        // fallback: treat as individual operators
        auto V_out = std::make_shared<QMPotential>(*this);
        mrcpp::deep_copy(*V_out, *this);
        out.push_back(O);
        out.push_back(V_out);
    }
    return out;
}

/** @brief compute real part of output
 *
 * @param inp: input orbital
 * @param dagger: apply complex conjugate potential
 *
 * Multiplies the potential with the orbital. The initial output grid is a
 * copy of the input orbital grid but NOT a copy of the potential grid.
 */
void QMPotential::calc(mrcpp::CompFunction<3> &out, mrcpp::CompFunction<3> &inp, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;
    if (out.Ncomp() > 0) MSG_ABORT("Output not empty");
    if (out.isShared()) MSG_ABORT("Cannot share this function");
    if (dagger and inp.iscomplex()) MSG_ERROR("Not implemented");
    if (inp.conjugate()) MSG_ERROR("Not implemented");

    mrcpp::CompFunction<3> &V = *this;
    double coef = 1.0;
    mrcpp::copy_grid(out, inp);
    mrcpp::multiply(prec, out, coef, inp, V, adap);
}

} // namespace mrchem
