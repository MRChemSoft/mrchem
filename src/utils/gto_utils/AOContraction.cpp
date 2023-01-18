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

using mrcpp::GaussExp;
using mrcpp::GaussFunc;

namespace mrchem {
namespace gto_utils {

// clang-format off
static const int s_gto[][3] = {
    {0, 0, 0}
};

int p_gto[][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
};

static const int d_gto[][3] = {
    {2, 0, 0},
    {1, 1, 0},
    {1, 0, 1},
    {0, 2, 0},
    {0, 1, 1},
    {0, 0, 2}
};

static const int f_gto[][3] = {
    {3, 0, 0},
    {2, 1, 0},
    {2, 0, 1},
    {1, 2, 0},
    {1, 1, 1},
    {1, 0, 2},
    {0, 3, 0},
    {0, 2, 1},
    {0, 1, 2},
    {0, 0, 3}
};

struct gtoDef {
    int l;
    int ncomp;
    const int *comp[(MAX_L + 1) * (MAX_L + 2) / 2];

};

static const gtoDef GTOS[MAX_L+1] = {
    { 0,  1, {s_gto[0]} },
    { 1,  3, {p_gto[0], p_gto[1], p_gto[2]}},
    { 2,  6, {d_gto[0], d_gto[1], d_gto[2], d_gto[3], d_gto[4], d_gto[5]}},
    { 3, 10, {f_gto[0], f_gto[1], f_gto[2], f_gto[3], f_gto[4], f_gto[5], f_gto[6], f_gto[7], f_gto[8], f_gto[9]}}
};
// clang-format on

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
GaussExp<3> AOContraction::getContraction(int m, const mrcpp::Coord<3> &center) const {
    assert(m >= 0 and m < this->nComp);
    GaussExp<3> ctr;
    double normFac = 1.0;
    const int *angMom = GTOS[this->L].comp[m];
    for (int i = 0; i < 3; i++) {
        switch (angMom[i]) {
            case 0:
                normFac *= 1.0;
                break;
            case 1:
                normFac *= 1.0;
                break;
            case 2:
                normFac *= 3.0;
                break;
            case 3:
                normFac *= 15.0;
                break;
            default:
                MSG_ERROR("We don't support g-functions at the moment");
        }
    }
    normFac = std::sqrt(normFac);

    std::array<int, 3> pow{angMom[0], angMom[1], angMom[2]};
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
