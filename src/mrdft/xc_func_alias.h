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

#include <vector>
#include <string>
#include <algorithm>
namespace mrdft {

/**
 * @brief Maps a functional name string (e.g., "PBE0", "LDA" or "XC_LDA_X", XC_GGA_X_B88) 
 * to its corresponding Libxc IDs and scaling coefficients
 * @note The input `name` is transformed to uppercase internally, making the
 * search case-insensitive
 * @param[in] name    Name of the functional
 * @param[in] ids     Vector to be populated with the IDs used by Libxc
 * @param[in] coefs   Vector to be populated with the corresponding scaling coefficients
 * @throws MSG_ABORT If the name is not a recognized internal shorthand and
 * is not found within the Libxc library
 * @example
 * std::vector<int> ids;
 * std::vector<double> coefs;
 * MapFuncName("LDA", ids, coefs);
 * // ids: {XC_LDA_C_VWN, XC_LDA_X}, coefs: {1.0, 1.0}
 */
void mapFunctionalName(std::string name, std::vector<int> &ids, std::vector<double> &coefs);

} // namespace mrdft
