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

#include <memory>

#include <MRCPP/Printer>

#include "mrchem.h"
#include "qmfunctions/Orbital.h"
#include "tensor/tensor_fwd.h"

/** @class QMOperator
 *
 * @brief Fundamental quantum mechanical operator
 *
 * Base class to handle operators and their application in the QM sense. Used to
 * build more complicated operators through the TensorOperator classes. This class
 * hierarchy should NOT be used directly, as the most important functionality is
 * protected. A proper interface is provided through RankZeroOperator.
 *
 * Note that operators are treated as real, with possibly a Complex scalar factor,
 * this factor is applied only if the function the ooperator is applied to is complex,
 * if the function is real, only the soft factor (out.func_ptr->data.c1[0]) is multiplied.
 *
 * Notes on naming conventions of derived operator classes:
 * Direct decendants of QMOperator should START with "QM", like QMPotential, QMSpin,
 * QMMomentum. Further decendants of QMPotential should END with "Potential", like
 * PositionPotential, NuclearPotential, CoulombPotential. Decendants of the
 * TensorOperators should end with "Operator" (except the perturbation operators
 * H_E_dip, H_B_dip, etc). E.i. the NuclearOperator IS a TensorOperator that
 * CONTAINS a NuclearPotential which IS a QMPotential which IS a QMOperator.
 * Capiche?
 *
 */
namespace mrchem {

class QMOperator {
public:
    QMOperator() = default;
    virtual ~QMOperator() {
        if (prec() > 0.0) MSG_ERROR("QMOperator not properly cleared");
    }

    double prec() { return this->apply_prec; }
    bool isImag() { return this->imag; }

    friend RankZeroOperator;
    virtual void setup(double prec) { setApplyPrec(prec); }

    protected:
    double apply_prec{-1.0};

    void setApplyPrec(double prec) {
        std::cout << "QMOperator::setApplyPrec -- Setting apply precision to " << prec << std::endl;
        if (this->apply_prec < 0.0) {
            std::cout << "QMOperator::setApplyPrec -- tut" << std::endl;
            this->apply_prec = prec;
            std::cout << "QMOperator::setApplyPrec -- pouet" << std::endl;
        } else if (not isSetup(prec)) {
            MSG_ERROR("Clear operator before setup with different prec!");
        }
    }
    void clearApplyPrec() { this->apply_prec = -1.0; }

    bool isSetup(double prec) const {
        MSG_INFO("QMOperator::setApplyPrec -- start ");
        double dPrec = std::abs(this->apply_prec - prec);
        double thrs = mrcpp::MachineZero;
        // MSG_INFO("QMOperator::setApplyPrec -- end ", " apply_prec=", this->apply_prec, " prec=", prec, " dPrec=", dPrec, " thrs=", thrs);
        std::cout << "QMOperator::setApplyPrec -- end " << " apply_prec=" << this->apply_prec << " prec=" << prec << " dPrec=" << dPrec << " thrs=" << thrs << std::endl;
        return (dPrec < thrs) ? true : false;
    }

    virtual void clear() { clearApplyPrec(); }

    virtual ComplexDouble evalf(const mrcpp::Coord<3> &r) const = 0;

    virtual Orbital apply(Orbital inp) = 0;
    virtual Orbital dagger(Orbital inp) = 0;
    virtual QMOperatorVector apply(std::shared_ptr<QMOperator> &O) = 0; //QMOperatorVector is an alias defined in tensor_fwd.h
    bool imag = false; // add imaginary unit prefactor, for faster application
};

} // namespace mrchem
