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

#include <vector>

#include "MRCPP/Gaussians"

#include "mrchem.h"

namespace mrchem {

namespace gto_utils {
class Intgrl;

class OrbitalExp final {
public:
    OrbitalExp(Intgrl &intgrl);
    ~OrbitalExp();

    int size() const { return this->orbitals.size(); }
    int getAngularMomentum(int n) const;

    mrcpp::GaussExp<3> getAO(int i) const { return *this->orbitals[i]; }
    mrcpp::GaussExp<3> getMO(int i, const DoubleMatrix &M) const;
    mrcpp::GaussExp<3> getDens(const DoubleMatrix &D) const;

    void rotate(const DoubleMatrix &U);

protected:
    bool cartesian;
    std::vector<mrcpp::GaussExp<3> *> orbitals;

    void readAOExpansion(Intgrl &intgrl);
    void transformToSpherical();
};

} // namespace gto_utils
} // namespace mrchem
