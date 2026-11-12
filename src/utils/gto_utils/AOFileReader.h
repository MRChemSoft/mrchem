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

#include "chemistry/Nucleus.h"
#include "utils/gto_utils/AOBasis.h"

#include <string>
#include <memory>

namespace mrchem {
namespace gto_utils {

/**
 * @brief Abstract base class for reading atomic orbital basis sets from files.
 *
 * This class provides an interface for reading atomic orbital basis sets from files.
 * Derived classes should implement the readBasisFile method to handle specific file formats.
 */
class AOFileReader {
public:
    virtual ~AOFileReader() = default;

    /**
     * @brief Get the number of nuclei in the basis set.
     * @return The number of nuclei.
     */
    size_t getNNuclei() const { return this->nuclei.size(); }

    /**
     * @brief Get the Nucleus object for the specified index.
     * @param i The index of the nucleus.
     * @return A shared pointer to the Nucleus object.
     */
    std::shared_ptr<Nucleus> getNucleus(size_t i) const { return this->nuclei[i]; }

    /**
     * @brief Get the AOBasis object for the specified index.
     * @param i The index of the basis set.
     * @return A shared pointer to the AOBasis object.
     */
    std::shared_ptr<AOBasis> getAOBasis(size_t i) const { return this->basis[i]; }
    
    /**
     * @brief Get the molecular orbital coefficients from a file.
     * @param mo_file The path to the molecular orbital file.
     * @return A DoubleMatrix containing the molecular orbital coefficients.
     */
    virtual DoubleMatrix getMOCoefficients(const std::string &MOFile) const = 0;
    
protected:
    /**
     * @brief Read the basis set.
     *
     * This method should be implemented by derived classes to read the basis set
     * from the specified file and populate the aoBasis member variable.
     */
    virtual void readBasisFile() = 0;

    std::vector<std::shared_ptr<Nucleus>> nuclei;           ///< Nucleus objects representing the nuclei in the basis set
    std::vector<std::shared_ptr<AOBasis>> basis; ///< AOBasis objects representing the basis sets for each atom type
    std::string AOFile;                                     ///< Path to the basis file
};

} // namespace gto_utils
} // namespace mrchem