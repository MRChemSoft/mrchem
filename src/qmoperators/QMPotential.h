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

#include "MRCPP/MWFunctions"

#include "QMOperator.h"

/** @class QMPotential
 *
 * @brief Operator defining a multiplicative potential
 *
 * Inherits the general features of a complex function from mrcpp::ComplexFunction and
 * implements the multiplication of this function with an Orbital. The actual
 * function representing the operator needs to be implemented in the derived
 * classes, where the *re and *im FunctionTree pointers should be assigned in
 * the setup() function and deallocated in the clear() function.
 *
 */

namespace mrchem {

class QMPotential : public mrcpp::ComplexFunction, public QMOperator {
public:
    explicit QMPotential(int adap, bool shared = false);
    QMPotential(const QMPotential &pot);
    QMPotential &operator=(const QMPotential &pot) = delete;
    virtual ~QMPotential() = default;

protected:
    int adap_build;

    ComplexDouble evalf(const mrcpp::Coord<3> &r) const override {
        ComplexDouble out(0.0, 0.0), i(0.0, 1.0);
        if (this->hasReal()) out += this->real().evalf(r);
        if (this->hasImag()) out += i * this->imag().evalf(r);
        return out;
    }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
    QMOperatorVector apply(std::shared_ptr<QMOperator> &O) override;

    void calcRealPart(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction &inp, bool dagger);
    void calcImagPart(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction &inp, bool dagger);
};

} // namespace mrchem
