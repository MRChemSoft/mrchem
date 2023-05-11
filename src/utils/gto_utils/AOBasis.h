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

#include "AOContraction.h"
#include <vector>

namespace mrchem {
namespace gto_utils {

class AOBasis final {
public:
    AOBasis();
    AOBasis(const AOBasis &bas);
    AOBasis &operator=(const AOBasis &bas);
    ~AOBasis();

    void append(const AOContraction &ctr);

    mrcpp::GaussExp<3> getAO(int n, const mrcpp::Coord<3> &center) const;
    mrcpp::GaussExp<3> getBasis(const mrcpp::Coord<3> &center) const;
    mrcpp::GaussExp<3> getNormBasis(const mrcpp::Coord<3> &center) const;

    AOContraction &getContraction(int n) { return *this->ctrs[n]; }
    const AOContraction &getContraction(int n) const { return *this->ctrs[n]; }

    int size() const { return this->ctrs.size(); }
    int getNFunc() const { return this->nFunc; }

    // This should print shell by shell
    friend std::ostream &operator<<(std::ostream &o, const AOBasis &b) {
        o << "    nFunc " << b.nFunc << std::endl;
        for (auto ctr : b.ctrs) o << "    " << *ctr;
        return o;
    }

private:
    int nPrim; ///< Total number of primitives in set
    int nFunc; ///< Total number of functions (all l-components included)
    std::vector<AOContraction *> ctrs;
};

} // namespace gto_utils
} // namespace mrchem
