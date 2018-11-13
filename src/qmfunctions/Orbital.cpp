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

#include "Orbital.h"
#include "orbital_utils.h"
#include "qmfunction_utils.h"

namespace mrchem {

/** @brief Default constructor
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 */
Orbital::Orbital()
        : QMFunction(false, nullptr, nullptr)
        , orb_data({-1, 0, 0, 1.0}) {}

/** @brief Constructor
 *
 * @param spin: electron spin (SPIN::Alpha/Beta/Paired)
 * @param occ: occupancy
 * @param rank: MPI ownership (-1 means all MPI ranks)
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 */
Orbital::Orbital(int spin, int occ, int rank)
        : QMFunction(false, nullptr, nullptr)
        , orb_data({rank, spin, occ, 1.0}) {
    if (this->spin() < 0) INVALID_ARG_ABORT;
    if (this->occ() < 0) {
        if (this->spin() == SPIN::Paired) this->orb_data.occ = 2;
        if (this->spin() == SPIN::Alpha) this->orb_data.occ = 1;
        if (this->spin() == SPIN::Beta) this->orb_data.occ = 1;
    }
}

/** @brief Copy constructor
 *
 * @param orb: orbital to copy
 *
 * Shallow copy: meta data is copied along with the *re and *im pointers,
 * NO transfer of ownership.
 */
Orbital::Orbital(const Orbital &orb)
        : QMFunction(orb)
        , orb_data(orb.orb_data) {}

/** @brief Assignment operator
 *
 * @param orb: orbital to copy
 *
 * Shallow copy: meta data is copied along with the *re and *im pointers,
 * NO transfer of ownership.
 */
Orbital &Orbital::operator=(const Orbital &orb) {
    if (this != &orb) {
        if (orb.isShared()) MSG_FATAL("Cannot shallow copy shared trees");
        this->orb_data = orb.orb_data;
        this->func_data = orb.func_data;
        this->re = orb.re;
        this->im = orb.im;
    }
    return *this;
}

/** @brief Parameter copy
 *
 * Returns a new orbital with the same spin, occupancy and rank_id as *this orbital.
 */
Orbital Orbital::paramCopy() const {
    return Orbital(this->spin(), this->occ(), this->rankID());
}

/** @brief Deep copy
 *
 * Returns a new orbital which is a full blueprint copy of *this orbital. This is
 * achieved by building a new grid for the real and imaginary parts and adding
 * in place.
 */
Orbital Orbital::deepCopy() {
    Orbital out(*this); // Shallow copy (should copy all meta data)
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
    return out; // Return shallow copy
}

/** @brief Complex conjugation
 *
 * Returns a new orbital which is a shallow copy of *this orbital, with a flipped
 * conjugate parameter. Pointer ownership is not transferred, so *this and the output
 * orbital points to the same MW representations of the real and imaginary parts,
 * however, they interpret the imaginary part with opposite sign.
 */
Orbital Orbital::dagger() const {
    Orbital out(*this); // Shallow copy
    out.func_data.conjugate = not this->func_data.conjugate;
    return out; // Return shallow copy
}

/** @brief Write orbital to disk
 *
 * @param file: file name prefix
 *
 * Given a file name prefix (e.g. "phi_0"), this will produce separate
 * binary files for meta data ("phi_0.meta"), real ("phi_0_re.tree")
 * and imaginary ("phi_0_im.tree") parts.
 */
void Orbital::saveOrbital(const std::string &file) {
    //writing meta data
    std::stringstream metafile;
    metafile << file << ".meta";

    //this flushes tree sizes
    FunctionData &func_data = getFunctionData();
    OrbitalData &orb_data = getOrbitalData();

    std::fstream f;
    f.open(metafile.str(), std::ios::out | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");
    f.write((char *)&func_data, sizeof(FunctionData));
    f.write((char *)&orb_data, sizeof(OrbitalData));
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

/** @brief Read orbital from disk
 *
 * @param file: file name prefix
 *
 * Given a file name prefix (e.g. "phi_0"), this will read separate
 * binary files for meta data ("phi_0.meta"), real ("phi_0_re.tree")
 * and imaginary ("phi_0_im.tree") parts.
 */
void Orbital::loadOrbital(const std::string &file) {
    if (hasReal()) MSG_ERROR("Orbital not empty");
    if (hasImag()) MSG_ERROR("Orbital not empty");

    //reading meta data
    std::stringstream fmeta;
    fmeta << file << ".meta";

    //this flushes tree sizes
    FunctionData &func_data = getFunctionData();
    OrbitalData &orb_data = getOrbitalData();

    std::fstream f;
    f.open(fmeta.str(), std::ios::in | std::ios::binary);
    if (f.is_open()) f.read((char *)&func_data, sizeof(FunctionData));
    if (f.is_open()) f.read((char *)&orb_data, sizeof(OrbitalData));
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

/** @brief Returns a character representing the spin (a/b/p) */
char Orbital::printSpin() const {
    char sp = 'u';
    if (this->spin() == SPIN::Paired) sp = 'p';
    if (this->spin() == SPIN::Alpha) sp = 'a';
    if (this->spin() == SPIN::Beta) sp = 'b';
    return sp;
}

/** @brief Pretty output of orbital meta data */
std::ostream &Orbital::print(std::ostream &o) const {
    int oldprec = mrcpp::Printer::setPrecision(12);
    o << std::setw(6) << this->rankID();
    o << std::setw(25) << this->norm();
    o << std::setw(5) << this->printSpin();
    o << std::setw(4) << this->occ();
    mrcpp::Printer::setPrecision(5);
    o << std::setw(15) << this->error();
    mrcpp::Printer::setPrecision(oldprec);
    return o;
}

} //namespace mrchem
