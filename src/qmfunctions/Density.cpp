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

#include <fstream>

#include "MRCPP/Printer"

#include "Density.h"
#include "mrchem.h"

namespace mrchem {

/** @brief Assignment operator
 *
 * @param dens: density to copy
 *
 * Shallow copy: meta data is copied along with the *re and *im pointers,
 * NO transfer of ownership.
 */
Density &Density::operator=(const Density &dens) {
    if (this != &dens) mrcpp::ComplexFunction::operator=(dens);
    return *this;
}

/** @brief Write density to disk
 *
 * @param file: file name prefix
 *
 * Given a file name prefix (e.g. "phi_0"), this will produce separate
 * binary files for meta data ("phi_0.meta"), real ("phi_0_re.tree")
 * and imaginary ("phi_0_im.tree") parts.
 */
void Density::saveDensity(const std::string &file) {
    // writing meta data
    std::stringstream metafile;
    metafile << file << ".meta";

    // this flushes tree sizes
    mrcpp::FunctionData &func_data = mrcpp::ComplexFunction::getFunctionData();

    std::fstream f;
    f.open(metafile.str(), std::ios::out | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");
    f.write((char *)&func_data, sizeof(mrcpp::FunctionData));
    f.close();

    // writing real part
    if (hasReal()) {
        std::stringstream fname;
        fname << file << "_re";
        real().saveTree(fname.str());
    }

    // writing imaginary part
    if (hasImag()) {
        std::stringstream fname;
        fname << file << "_im";
        imag().saveTree(fname.str());
    }
}

/** @brief Read density from disk
 *
 * @param file: file name prefix
 *
 * Given a file name prefix (e.g. "phi_0"), this will read separate
 * binary files for meta data ("phi_0.meta"), real ("phi_0_re.tree")
 * and imaginary ("phi_0_im.tree") parts.
 */
void Density::loadDensity(const std::string &file) {
    if (hasReal()) MSG_ERROR("Density not empty");
    if (hasImag()) MSG_ERROR("Density not empty");

    // reading meta data
    std::stringstream fmeta;
    fmeta << file << ".meta";

    // this flushes tree sizes
    mrcpp::FunctionData &func_data = mrcpp::ComplexFunction::getFunctionData();

    std::fstream f;
    f.open(fmeta.str(), std::ios::in | std::ios::binary);
    if (f.is_open()) f.read((char *)&func_data, sizeof(mrcpp::FunctionData));
    f.close();

    // reading real part
    if (func_data.real_size > 0) {
        std::stringstream fname;
        fname << file << "_re";
        alloc(NUMBER::Real);
        real().loadTree(fname.str());
    }

    // reading imaginary part
    if (func_data.imag_size > 0) {
        std::stringstream fname;
        fname << file << "_im";
        alloc(NUMBER::Imag);
        imag().loadTree(fname.str());
    }
}

} // namespace mrchem
