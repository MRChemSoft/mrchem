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
#include <vector>
#include <H5Cpp.h>

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
     * @param nMOs The number of molecular orbitals to read.
     * @param nSpinors The number of spinors to read.
     * @return A vector of ComplexMatrix containing the spinor molecular orbital coefficients.
     */
    std::vector<ComplexMatrix> getMOCoefficients(const std::string &MOFile, size_t nMOs, size_t nSpinors) const override final;

    /**
     * @brief Get the small component AOBasis object for the specified index.
     * @param i The index of the basis set.
     * @return A shared pointer to the small component AOBasis object, or nullptr if not available.
     */
    std::shared_ptr<AOBasis> getSmallAOBasis(size_t i) const override final;

private:
    /**
     * @brief Read the basis set in Dirac format.
     */
    void readBasisFile() override final;

    /**
     * @brief Read the atoms in the basis set.
     * @param file The HDF5 file object.
     */
    void readAtomBlock(const H5::H5File &file);

    /**
     * @brief Read one DIRAC AO basis group (/input/aobasis/1 = large, /input/aobasis/2 = small).
     * @param g The HDF5 group of the basis.
     * @return One AOBasis per nucleus, in the order of this->nuclei.
     */
    std::vector<std::shared_ptr<AOBasis>> readAOBasis(const H5::Group &g) const;

    /**
     * @brief Find the nucleus at a given position.
     * @param r Position in bohr.
     * @return Index into this->nuclei.
     */
    size_t findAtom(const mrcpp::Coord<3> &r) const;

    std::vector<std::shared_ptr<AOBasis>> smallBasis; ///< AOBasis objects representing the small component basis sets for each atom type
};

} // namespace gto_utils
} // namespace mrchem