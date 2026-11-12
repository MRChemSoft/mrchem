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

#include <string>
#include <memory>

namespace mrchem {
class OrbitalVector;
namespace gto_utils {
        
/**
 * @brief Abstract base class for reading atomic orbital basis sets from files.
 *
 * This class provides an interface for reading atomic orbital basis sets from files.
 * Derived classes should implement the readBasisFile method to handle specific file formats.
 */
class GTOFileReader {
public:
    /**
     * @brief Construct a new GTOFileReader object.
     * @param GTOType The type of GTOs to use (e.g., "Dalton").
     * @param basisFile The path to the atomic orbital basis set file.
     */
    GTOFileReader(const std::string &GTOType, const std::string &basisFile);
    virtual ~GTOFileReader() = default;

    /**
     * @brief Project molecular orbitals from a file.
     * @param Phi The vector of molecular orbitals.
     * @param mo_file The path to the molecular orbital file.
     * @param prec The precision for the projection.
     * @param screen The screening parameter for the GTOs.
     */
    void projectMO(OrbitalVector &Phi, const std::string &mo_file, double prec, double screen);
private:

    std::unique_ptr<AOFileReader> reader;
};

} // namespace gto_utils
} // namespace mrchem