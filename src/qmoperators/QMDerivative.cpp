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

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "QMDerivative.h"

#include "QMIdentity.h"
#include "QMPotential.h"
#include "QMSpin.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

using QMOperator_p = std::shared_ptr<mrchem::QMOperator>;

namespace mrchem {

QMDerivative::QMDerivative(int d, std::shared_ptr<DerivativeOperator<3>> D, bool im)
        : QMOperator()
        , imag(im)
        , apply_dir(d)
        , derivative(D) {}

QMDerivative::QMDerivative(const QMDerivative &inp)
        : QMOperator()
        , imag(inp.imag)
        , apply_dir(inp.apply_dir)
        , derivative(inp.derivative) {}

Orbital QMDerivative::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == nullptr) MSG_ERROR("No derivative operator");

    auto dir = this->apply_dir;
    auto &D = *this->derivative;

    Orbital out = inp.paramCopy();
    if (this->isReal()) {
        if (inp.hasReal()) {
            out.alloc(NUMBER::Real);
            mrcpp::apply(out.real(), D, inp.real(), dir);
        }
        if (inp.hasImag()) {
            out.alloc(NUMBER::Imag);
            mrcpp::apply(out.imag(), D, inp.imag(), dir);
            if (inp.conjugate()) out.imag().rescale(-1.0);
        }
    } else {
        if (inp.hasImag()) {
            out.alloc(NUMBER::Real);
            mrcpp::apply(out.real(), D, inp.imag(), dir);
            if (inp.conjugate()) out.real().rescale(-1.0);
        }
        if (inp.hasReal()) {
            out.alloc(NUMBER::Imag);
            mrcpp::apply(out.imag(), D, inp.real(), dir);
            out.imag().rescale(-1.0);
        }
    }
    return out;
}

Orbital QMDerivative::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

QMOperatorVector QMDerivative::apply(QMOperator_p &O) {
    QMIdentity *I = dynamic_cast<QMIdentity *>(&(*O));
    QMPotential *V_inp = dynamic_cast<QMPotential *>(&(*O));

    QMOperatorVector out;
    if (I) {
        // O == identity: skip it
        out.push_back(std::make_shared<QMDerivative>(*this));
    } else if (V_inp) {
        // O == potential: compute its derivative
        auto &D = *this->derivative;
        auto d = this->apply_dir;
        auto V_out = std::make_shared<QMPotential>(*V_inp);
        if (this->isReal()) {
            if (V_inp->hasReal()) {
                V_out->alloc(NUMBER::Real);
                mrcpp::apply(V_out->real(), D, V_inp->real(), d);
            }
            if (V_inp->hasImag()) {
                V_out->alloc(NUMBER::Imag);
                mrcpp::apply(V_out->imag(), D, V_inp->imag(), d);
                if (V_inp->conjugate()) V_out->imag().rescale(-1.0);
            }
        } else {
            if (V_inp->hasImag()) {
                V_out->alloc(NUMBER::Real);
                mrcpp::apply(V_out->real(), D, V_inp->imag(), d);
                if (V_inp->conjugate()) V_out->real().rescale(-1.0);
            }
            if (V_inp->hasReal()) {
                V_out->alloc(NUMBER::Imag);
                mrcpp::apply(V_out->imag(), D, V_inp->real(), d);
                V_out->imag().rescale(-1.0);
            }
        }
        out.push_back(V_out);
    } else {
        // fallback: treat as individual operators
        out.push_back(O);
        out.push_back(std::make_shared<QMDerivative>(*this));
    }
    return out;
}

} // namespace mrchem
