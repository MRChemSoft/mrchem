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

#include "QMOperator.h"

namespace mrchem {

/** @class QMSpin
 *
 * @brief Cartesian spin operator
 *
 * Defined by:
 *
 *  S_{x} alpha = (1/2)*beta
 *  S_{y} alpha = (i/2)*beta
 *  S_{z} alpha = (1/2)*alpha
 *
 *  S_{x} beta =  (1/2)*alpha
 *  S_{y} beta = -(i/2)*alpha
 *  S_{z} beta = -(1/2)*beta
 *
 */

class QMSpin final : public QMOperator {
public:
    QMSpin(int d)
            : D(d) {}
    QMSpin(const QMSpin &S)
            : D(S.D) {}

private:
    const int D;

    ComplexDouble evalf(const mrcpp::Coord<3> &r) const override { return 0.0; }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
    QMOperatorVector apply(std::shared_ptr<QMOperator> &O) override;
};

/** @class QMAlpha
 *
 * @brief Spin alpha operator
 *
 * Operator that acts as Identity when applied on alpha (or paired) orbitals,
 * and as Zero when applid on beta orbitals.
 *
 * Can be defined as:
 *
 *  S_{alpha} = S_{+} * S_{-},
 *
 * using the shift operators:
 *
 *  S_{+} = S_{x} + i*S_{y}
 *  S_{-} = S_{x} - i*S_{y}
 *
 * Could have been implemented using the Cartesian spin operators above,
 * but we rather do it explicitly as Identity or Zero, to avoid unnecessary
 * deep copies and addition/subtraction of identical contributions.
 *
 */

class QMAlpha final : public QMOperator {
public:
    QMAlpha() = default;

private:
    ComplexDouble evalf(const mrcpp::Coord<3> &r) const override { return 0.0; }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
    QMOperatorVector apply(std::shared_ptr<QMOperator> &O) override;
};

/** @class QMBeta
 *
 * @brief Spin beta operator
 *
 * Operator that acts as Identity when applied on beta (or paired) orbitals,
 * and as Zero when applid on alpha orbitals.
 *
 * Can be defined as:
 *
 *  S_{beta} = S_{-} * S_{+},
 *
 * using the shift operators:
 *
 *  S_{+} = S_{x} + i*S_{y}
 *  S_{-} = S_{x} - i*S_{y}
 *
 * Could have been implemented using the Cartesian spin operators above,
 * but we rather do it explicitly as Identity or Zero, to avoid unnecessary
 * deep copies and addition/subtraction of identical contributions.
 *
 */

class QMBeta final : public QMOperator {
public:
    QMBeta() = default;

private:
    ComplexDouble evalf(const mrcpp::Coord<3> &r) const override { return 0.0; }

    Orbital apply(Orbital inp) override;
    Orbital dagger(Orbital inp) override;
    QMOperatorVector apply(std::shared_ptr<QMOperator> &O) override;
};

} // namespace mrchem
