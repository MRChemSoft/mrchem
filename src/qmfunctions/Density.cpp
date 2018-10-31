/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

namespace mrchem {

/** @brief Assignment operator
 *
 * @param dens: density to copy
 *
 * Shallow copy: meta data is copied along with the *re and *im pointers,
 * NO transfer of ownership.
 */
Density& Density::operator=(const Density &dens) {
    if (this != &dens) {
        if (dens.isShared()) MSG_FATAL("Cannot shallow copy shared trees");
        this->func_data = dens.func_data;
        this->re = dens.re;
        this->im = dens.im;
    }
    return *this;
}

/** @brief Deep copy
 *
 * Returns a new density which is a full blueprint copy of *this density. This is
 * achieved by building a new grid for the real and imaginary parts and adding
 * in place.
 */
Density Density::deepCopy() {
    Density out(*this); // Shallow copy (should copy all meta data)
    out.clear();        // Remove *re and *im pointers
    if (this->hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::copy_grid(out.real(), this->real());
        mrcpp::copy_func(out.real(), this->real());
    }
    if (this->hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::copy_grid(out.imag(), this->imag());
        mrcpp::copy_func(out.imag(), this->imag());
    }
    return out;         // Return shallow copy
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
    //writing meta data
    std::stringstream metafile;
    metafile << file << ".meta";

    //this flushes tree sizes
    FunctionData &func_data = getFunctionData();

    std::fstream f;
    f.open(metafile.str(), std::ios::out | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");
    f.write((char *) &func_data, sizeof(FunctionData));
    f.close();

    //writing real part
    if (hasReal()) {
        std::stringstream fname;
        fname << file << "_re";
        this->real().saveTree(fname.str());
    }

    //writing imaginary part
    if (hasImag()) {
        std::stringstream fname;
        fname << file << "_im";
        this->imag().saveTree(fname.str());
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

    //reading meta data
    std::stringstream fmeta;
    fmeta << file << ".meta";

    //this flushes tree sizes
    FunctionData &func_data = getFunctionData();

    std::fstream f;
    f.open(fmeta.str(), std::ios::in | std::ios::binary);
    if (f.is_open()) f.read((char *) &func_data, sizeof(FunctionData));
    f.close();

    //reading real part
    if (func_data.nChunksReal > 0) {
        std::stringstream fname;
        fname << file << "_re";
        alloc(NUMBER::Real);
        this->real().loadTree(fname.str());
    }

    //reading imaginary part
    if (func_data.nChunksImag > 0) {
        std::stringstream fname;
        fname << file << "_im";
        alloc(NUMBER::Imag);
        this->imag().loadTree(fname.str());
    }
}

} //namespace mrchem

