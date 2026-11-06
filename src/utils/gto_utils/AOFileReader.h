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

#include "utils/gto_utils/AOBasis.h"

#include <string>
#include <memory>

namespace mrchem {

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
     * @brief Get the AOBasis object. If the basis file has not been read yet, this will read it.
     * @return shared pointer to AOBasis object
    */
    std::shared_ptr<gto_utils::AOBasis> getAOBasis();
    
protected:
    /**
     * @brief Read the basis set.
     *
     * This method should be implemented by derived classes to read the basis set
     * from the specified file and populate the aoBasis member variable.
     */
    virtual void readBasisFile() = 0;

    std::shared_ptr<gto_utils::AOBasis> aoBasis; ///< Shared pointer to the AOBasis object
    std::string AOFile; ///< Path to the basis file
};
} // namespace mrchem