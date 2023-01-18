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

#include <iostream>
#include <vector>

#include "MRCPP/Gaussians"

namespace mrchem {
namespace gto_utils {

static const int MAX_L = 3; // Max f-functions

/** Contracted atomic orbital.
 *
 * AOContraction defines a singel, contracted atomic orbital as a GaussExp
 * with coefficients.
 *
 */
class AOContraction final {
public:
    AOContraction(int l = 0);

    void append(double e, double c);
    mrcpp::GaussExp<3> getNormContraction(int m, const mrcpp::Coord<3> &center) const;
    mrcpp::GaussExp<3> getContraction(int m, const mrcpp::Coord<3> &center) const;

    int getNComp() const { return this->nComp; }
    int getMoment() const { return this->L; }
    void setMoment(int l) { this->L = l; }

    friend std::ostream &operator<<(std::ostream &o, const AOContraction &c) {
        o << " " << c.L << " " << c.coefs.size() << std::endl;
        for (unsigned int i = 0; i < c.expo.size(); i++) { o << "    " << c.expo[i] << "   " << c.coefs[i] << std::endl; }
        return o;
    }

protected:
    int L;
    int nComp;
    std::vector<double> expo;
    std::vector<double> coefs;
};

} // namespace gto_utils
} // namespace mrchem
