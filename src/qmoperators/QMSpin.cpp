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

#include "QMSpin.h"

#include "QMDerivative.h"
#include "QMIdentity.h"
#include "QMPotential.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using QMOperator_p = std::shared_ptr<mrchem::QMOperator>;

namespace mrchem {

///////////////////////////////////////////////
//////////          QMSpin           //////////
///////////////////////////////////////////////

Orbital QMSpin::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    ComplexDouble coef(0.0, 0.0);
    switch (inp.spin()) {
        case SPIN::Alpha:
            if (this->D == 0) coef = ComplexDouble(0.5, 0.0);
            if (this->D == 1) coef = ComplexDouble(0.0, 0.5);
            if (this->D == 2) coef = ComplexDouble(0.5, 0.0);
            break;
        case SPIN::Beta:
            if (this->D == 0) coef = ComplexDouble(0.5, 0.0);
            if (this->D == 1) coef = ComplexDouble(0.0, -0.5);
            if (this->D == 2) coef = ComplexDouble(-0.5, 0.0);
            break;
        default:
            MSG_ABORT("Cannot apply spin operator on paired orbital");
    }

    Orbital out = inp.paramCopy();
    qmfunction::deep_copy(out, inp);
    out.rescale(coef);

    // Flip spin for s_x and s_y
    if (this->D == 0 or this->D == 1) {
        if (inp.spin() == SPIN::Alpha) out.setSpin(SPIN::Beta);
        if (inp.spin() == SPIN::Beta) out.setSpin(SPIN::Alpha);
    }

    return out;
}

Orbital QMSpin::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

QMOperatorVector QMSpin::apply(QMOperator_p &O) {
    QMOperatorVector out;

    QMIdentity *I = dynamic_cast<QMIdentity *>(&(*O));
    if (I) {
        // O == identity: skip it
        out.push_back(std::make_shared<QMSpin>(*this));
    } else {
        // fallback: treat as individual operators
        out.push_back(O);
        out.push_back(std::make_shared<QMSpin>(*this));
    }

    return out;
}

///////////////////////////////////////////////
//////////          QMAlpha          //////////
///////////////////////////////////////////////

/** @brief Identity on alpha orbital, zero on beta orbital */
Orbital QMAlpha::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out(SPIN::Alpha, 1, inp.rankID());
    if (inp.spin() == SPIN::Alpha || inp.spin() == SPIN::Paired) qmfunction::deep_copy(out, inp);
    return out;
}

/** @brief Identity on alpha orbital, zero on beta orbital */
Orbital QMAlpha::dagger(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out(SPIN::Alpha, 1, inp.rankID());
    if (inp.spin() == SPIN::Alpha || inp.spin() == SPIN::Paired) qmfunction::deep_copy(out, inp);
    return out;
}

QMOperatorVector QMAlpha::apply(QMOperator_p &O) {
    QMIdentity *I = dynamic_cast<QMIdentity *>(&(*O));
    QMAlpha *A = dynamic_cast<QMAlpha *>(&(*O));

    QMOperatorVector out;
    if (I or A) {
        // O == identity: skip it
        // O == alpha: alpha^2 == alpha -> skip it
        out.push_back(std::make_shared<QMAlpha>(*this));
    } else {
        // fallback: treat as individual operators
        out.push_back(O);
        out.push_back(std::make_shared<QMAlpha>(*this));
    }

    return out;
}

///////////////////////////////////////////////
//////////          QMBeta           //////////
///////////////////////////////////////////////

/** @brief Identity on beta orbital, zero on alpha orbital */
Orbital QMBeta::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out(SPIN::Beta, 1, inp.rankID());
    if (inp.spin() == SPIN::Beta || inp.spin() == SPIN::Paired) qmfunction::deep_copy(out, inp);
    return out;
}

/** @brief Identity on beta orbital, zero on alpha orbital */
Orbital QMBeta::dagger(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out(SPIN::Beta, 1, inp.rankID());
    if (inp.spin() == SPIN::Beta || inp.spin() == SPIN::Paired) qmfunction::deep_copy(out, inp);
    return out;
}

QMOperatorVector QMBeta::apply(QMOperator_p &O) {
    QMOperatorVector out;

    QMIdentity *I = dynamic_cast<QMIdentity *>(&(*O));
    QMBeta *B = dynamic_cast<QMBeta *>(&(*O));
    if (I or B) {
        // O == identity: skip it
        // O == beta: beta^2 == beta -> skip it
        out.push_back(std::make_shared<QMBeta>(*this));
    } else {
        // fallback: treat as individual operators
        out.push_back(O);
        out.push_back(std::make_shared<QMBeta>(*this));
    }

    return out;
}

} // namespace mrchem
