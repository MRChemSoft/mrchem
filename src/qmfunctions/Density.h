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

#include "MRCPP/MWFunctions"

#include "Orbital.h"

/** @class Density
 *
 * @brief General complex-valued function to handle densities (incl trans. densities)
 *
 * Inherits the general features of a complex function from mrcpp::ComplexFunction which
 * means separate MW function representations for the real and imaginary parts.
 * Note that there are several options for copying/assignment: the proper copy
 * constructor and assignment operator are *shallow* copies, which means that
 * they simply copy the *re and *im pointers (no transfer of ownership).
 * Additionaly, there is a deepCopy() which returns a *full* copy of the density.
 *
 * NOTE: since the standard copies are shallow copies and several densitys can
 * point to the same MW functions, it is YOUR responibility to keep track of the
 * ownership and free the FunctionTree pointers before the last density goes out
 * of scope.
 */

namespace mrchem {

class Density final : public mrcpp::ComplexFunction {
public:
    explicit Density(bool share)
            : mrcpp::ComplexFunction(0, -1, -1, share) {}
    Density(const Density &dens)
            : mrcpp::ComplexFunction(dens) {}
    Density &operator=(const Density &dens);

    void saveDensity(const std::string &file);
    void loadDensity(const std::string &file);
};

} // namespace mrchem
