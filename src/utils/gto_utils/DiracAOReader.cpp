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

#include "utils/gto_utils/DiracAOReader.h"

#include <MRCPP/Printer>

#include <filesystem>
#include <H5Cpp.h>

namespace mrchem {
namespace gto_utils {

DiracAOReader::DiracAOReader(const std::string &basisFile) {
    if (not std::filesystem::is_regular_file(basisFile))
        MSG_ABORT("Failed to open basis set file: " << basisFile);
    this->AOFile = basisFile;
    readBasisFile();
}

void DiracAOReader::readBasisFile() {
    return;
}

DoubleMatrix DiracAOReader::getMOCoefficients(const std::string &MOFile) const {
    DoubleMatrix C;
    return C;
}

} // namespace gto_utils
} // namespace mrchem