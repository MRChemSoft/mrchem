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

class DaltonAOReader final : public AOFileReader {
public:
    /**
     * @brief Construct a new DaltonAOReader object.
     * @param file The path to the Dalton basis file.
     */
    DaltonAOReader(const std::string &basisFile);

private:
    /**
     * @brief Read the basis set in Dalton format.
     */
    void readBasisFile() override final;

    /**
     * @brief Read the basis set for a single atom type.
     * @param ifs The input file stream.
     */
    void readAtomBlock(std::ifstream &ifs);

    /**
     * @brief Read the atomic data for a single atom type.
     * @param ifs The input file stream.
     * @param nAtoms The number of atoms of this type.
     * @param z The atomic number of this atom type.
     */
    void readAtomData(std::ifstream &ifs, size_t nAtoms, double z);

    /**
     * @brief Read the contraction block for a single shell.
     * @param ifs The input file stream.
     * @param bas The AOBasis object to populate.
     * @param l The angular momentum quantum number of the shell.
     */
    void readContractionBlock(std::ifstream &ifs, AOBasis &bas, size_t l);
};

} // namespace gto_utils
} // namespace mrchem