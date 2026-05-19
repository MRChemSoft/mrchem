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

#include "AOContraction.h"

#include "utils/math_utils.h"

using mrcpp::GaussExp;
using mrcpp::GaussFunc;

namespace mrchem {
namespace gto_utils {

// TODO: Docs
static std::array<int, 3> ind_to_crt(int l, int i) {
    int lx = l;

    for (int j = 1; i > l - lx; j++) {
        lx--;
        i -= j;
    }

    int ly = l - lx - i;
    int lz = l - lx - ly;

    return {lx, ly, lz};
}

AOContraction::AOContraction(int l) {
    this->L = l;
    this->nComp = (this->L + 1) * (this->L + 2) / 2;
}

GaussExp<3> AOContraction::getNormContraction(int m, const mrcpp::Coord<3> &center) const {
    GaussExp<3> ctr = getContraction(m, center);
    ctr.normalize();
    return ctr;
}

/** Normalization goes like this (thanks Radovan)

    < AO | AO > =   1 for s, px, py, pz, dxy, dxz, dyz, fxyz
            3 for dxx, dyy, dzz, fxxy, ...
               15     fxxx, fyyy, fzzz, gxxxy, ...
              105     gxxxx, ...
              ...     ...
            9 for gxxyy, ...
              etc     ...
*/
double cartesianNormFac(int l) {
    return math_utils::double_factorial(2 * l - 1);
}

double cartesianNormFac(int lx, int ly, int lz) {
    int l = lx + ly + lz;

    double n = cartesianNormFac(lx) * cartesianNormFac(ly) * cartesianNormFac(lz) / cartesianNormFac(l);

    return std::sqrt(n);
}

GaussExp<3> AOContraction::getContraction(int m, const mrcpp::Coord<3> &center) const {
    assert(m >= 0 and m < this->nComp);
    GaussExp<3> ctr;
    std::array<int, 3> pow = ind_to_crt(this->L, m);
    double normFac = cartesianNormFac(pow[0], pow[1], pow[2]);

    for (unsigned int i = 0; i < expo.size(); i++) {
        GaussFunc<3> gto(this->expo[i], 1.0, center, pow);
        gto.normalize();
        gto.setCoef(normFac * gto.getCoef() * this->coefs[i]);
        ctr.append(gto);
    }
    return ctr;
}

void AOContraction::append(double e, double c) {
    this->expo.push_back(e);
    this->coefs.push_back(c);
}

} // namespace gto_utils
} // namespace mrchem
