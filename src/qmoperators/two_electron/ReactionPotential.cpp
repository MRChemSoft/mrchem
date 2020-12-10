/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "ReactionPotential.h"

#include "qmfunctions/qmfunction_utils.h"

using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

ReactionPotential::ReactionPotential(OrbitalVector_p Phi_p, SCRF helper)
        : QMPotential(1, false)
        , Phi(Phi_p)
        , helper(helper) {}

void ReactionPotential::setup(double prec) {
    setApplyPrec(prec);
    QMFunction &temp = *this;
    if (this->first_iteration) {
        this->first_iteration = false;
        return;
    }
    qmfunction::deep_copy(temp, this->helper.setup(prec, this->Phi));
}

void ReactionPotential::clear() {
    clearApplyPrec();
    this->helper.clear();
}

} // namespace mrchem
