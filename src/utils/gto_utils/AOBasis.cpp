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

#include "MRCPP/Printer"

#include "AOBasis.h"
#include "AOContraction.h"

using mrcpp::GaussExp;

namespace mrchem {
namespace gto_utils {

AOBasis::AOBasis() {
    this->nPrim = 0;
    this->nFunc = 0;
}

AOBasis::AOBasis(const AOBasis &bas) {
    this->nPrim = 0;
    this->nFunc = 0;
    for (int i = 0; i < bas.size(); i++) append(bas.getContraction(i));
}

AOBasis &AOBasis::operator=(const AOBasis &bas) {
    if (this != &bas) {
        this->nPrim = 0;
        this->nFunc = 0;
        for (int i = 0; i < bas.size(); i++) append(bas.getContraction(i));
    }
    return *this;
}

AOBasis::~AOBasis() {
    for (unsigned int i = 0; i < ctrs.size(); i++) {
        if (this->ctrs[i] != nullptr) delete this->ctrs[i];
    }
}

void AOBasis::append(const AOContraction &ctr) {
    this->ctrs.push_back(new AOContraction(ctr));
    this->nPrim++;
    this->nFunc += ctr.getNComp();
}

GaussExp<3> AOBasis::getAO(int n, const mrcpp::Coord<3> &center) const {
    assert(n >= 0 and n < nFunc);
    int m = 0;
    for (auto i : ctrs) {
        const AOContraction &ctr = *i;
        for (int j = 0; j < ctr.getNComp(); j++) {
            if (m == n) return ctr.getNormContraction(j, center);
            m++;
        }
    }
    MSG_ABORT("Something is terribly wrong");
}

GaussExp<3> AOBasis::getBasis(const mrcpp::Coord<3> &center) const {
    NOT_IMPLEMENTED_ABORT;
    GaussExp<3> abas;
    for (auto i : this->ctrs) {
        const AOContraction &ctr = *i;
        for (int m = 0; m < ctr.getNComp(); m++) abas.append(ctr.getContraction(m, center));
    }
    return abas;
}

GaussExp<3> AOBasis::getNormBasis(const mrcpp::Coord<3> &center) const {
    NOT_IMPLEMENTED_ABORT;
    GaussExp<3> abas;
    for (auto i : this->ctrs) {
        const AOContraction &ctr = *i;
        for (int m = 0; m < ctr.getNComp(); m++) abas.append(ctr.getNormContraction(m, center));
    }
    return abas;
}

} // namespace gto_utils
} // namespace mrchem
