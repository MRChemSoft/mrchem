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

#include "utils/gto_utils/AOFileReader.h"

#include <fstream>

namespace mrchem {
namespace gto_utils {

class DiracAOReader final : public AOFileReader {
public:
    /**
     * @brief Construct a new DiracAOReader object.
     * @param file The path to the Dirac basis file.
     */
    DiracAOReader(const std::string &basisFile);

    /**
     * @brief Get the molecular orbital coefficients from a file.
     * @param MOFile The path to the molecular orbital file.
     * @return A DoubleMatrix containing the molecular orbital coefficients.
     */
    DoubleMatrix getMOCoefficients(const std::string &MOFile) const override final;

private:
    /**
     * @brief Read the basis set in Dalton format.
     */
    void readBasisFile() override final;
};

} // namespace gto_utils
} // namespace mrchem