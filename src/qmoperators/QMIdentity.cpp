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

#include "QMIdentity.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using QMOperator_p = std::shared_ptr<mrchem::QMOperator>;

namespace mrchem {

/** Identity operator is a deep copy */
Orbital QMIdentity::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    Orbital out = inp.paramCopy();
    mrcpp::cplxfunc::deep_copy(out, inp);
    return out;
}

/** Identity operator is a deep copy */
Orbital QMIdentity::dagger(Orbital inp) {
    return apply(inp);
}

QMOperatorVector QMIdentity::apply(QMOperator_p &O) {
    // this == identity: always skip it and pass O through
    QMOperatorVector out;
    out.push_back(O);
    return out;
}

} // namespace mrchem
