/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "MRCPP/MWFunctions"

#include "chemistry/Nucleus.h"

namespace mrchem {

class NuclearGradientFunction : public mrcpp::RepresentableFunction<3> {
public:
    NuclearGradientFunction(int d, const Nucleus &nuc, double c)
        : dir(d), smooth(c), nucleus(nuc) { }

    double evalf(const double *r) const;

    Nucleus &getNucleus() { return this->nucleus; }
    const Nucleus &getNucleus() const { return this->nucleus; }

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

protected:
    int dir;
    double smooth;
    Nucleus nucleus;

    double du_dr(double r1) const;
};

} //namespace mrchem
