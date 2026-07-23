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

#include "utils/PolyInterpolator.h"
#include <Eigen/Dense>
#include <vector>

class HirshfeldRadInterpolater {

public:
    /**
     * @brief Construct a new Rad Interpolater object
     * @param element The element for which the ZORA potential is to be interpolated
     * @param data_dir The directory containing the ZORA potential data
     */
    HirshfeldRadInterpolater(const std::string element, std::string data_dir, bool writeToFile = false);

    /**
     * @brief Evaluate the interpolated function at a given point
     * @param r The point at which to evaluate
     * @return The interpolated value at the given point
     */
    double evalf(const double &r) const;

    /**
     * @brief Integrates the atomic charge density to get the charge of the tabulated
     * atomic density. The numeric is performed from 0 to 20 bohr. Only useful for debugging.
     */
    double getNorm() const;

protected:
    /**
     * @brief The interpolator for the atomic density
     */
    std::shared_ptr<interpolation_utils::PolyInterpolator> lnRho = nullptr;
    /**
     * @brief Write the interpolated density to a file for debugging.
     * @param path The path to the file to write the interpolated density to.
     */
    void writeInterpolatedDensity(const std::string path);
};
