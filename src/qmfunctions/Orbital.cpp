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

#include "Orbital.h"
#include "orbital_utils.h"
#include "qmfunction_utils.h"

namespace mrchem {

/** @brief Default constructor
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 */
Orbital::Orbital()
        : QMFunction(false)
        , orb_data({-1, 0, 0}) {}

/** @brief Constructor
 *
 * @param spin: electron spin (SPIN::Alpha/Beta/Paired)
 * @param occ: occupation
 * @param rank: MPI ownership (-1 means all MPI ranks)
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 */
Orbital::Orbital(int spin, int occ, int rank)
        : QMFunction(false)
        , orb_data({rank, spin, occ}) {
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
        QMFunction::operator=(orb);
        this->orb_data = orb.orb_data;
    }
    return *this;
}

Orbital &Orbital::operator=(const QMFunction &func) {
    if (this != &func) QMFunction::operator=(func);
    return *this;
}

/** @brief Parameter copy
 *
 * Returns a new orbital with the same spin, occupation and rank_id as *this orbital.
 */
Orbital Orbital::paramCopy() const {
    return Orbital(this->spin(), this->occ(), this->rankID());
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
    out.conj = not this->conjugate();
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
    // writing meta data
    std::stringstream metafile;
    metafile << file << ".meta";

    // this flushes tree sizes
    FunctionData &func_data = getFunctionData();
    OrbitalData &orb_data = getOrbitalData();

    std::fstream f;
    f.open(metafile.str(), std::ios::out | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");
    f.write((char *)&func_data, sizeof(FunctionData));
    f.write((char *)&orb_data, sizeof(OrbitalData));
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

    // reading meta data
    std::stringstream fmeta;
    fmeta << file << ".meta";

    // this flushes tree sizes
    FunctionData &func_data = getFunctionData();
    OrbitalData &orb_data = getOrbitalData();

    std::fstream f;
    f.open(fmeta.str(), std::ios::in | std::ios::binary);
    if (f.is_open()) f.read((char *)&func_data, sizeof(FunctionData));
    if (f.is_open()) f.read((char *)&orb_data, sizeof(OrbitalData));
    f.close();

    std::array<int, 3> corner{func_data.corner[0], func_data.corner[1], func_data.corner[2]};
    std::array<int, 3> boxes{func_data.boxes[0], func_data.boxes[1], func_data.boxes[2]};
    mrcpp::BoundingBox<3> world(func_data.scale, corner, boxes);

    mrcpp::MultiResolutionAnalysis<3> *mra = nullptr;
    if (func_data.type == mrcpp::Interpol) {
        mrcpp::InterpolatingBasis basis(func_data.order);
        mra = new mrcpp::MultiResolutionAnalysis<3>(world, basis, func_data.depth);
    } else if (func_data.type == mrcpp::Legendre) {
        mrcpp::LegendreBasis basis(func_data.order);
        mra = new mrcpp::MultiResolutionAnalysis<3>(world, basis, func_data.depth);
    } else {
        MSG_ABORT("Invalid basis type!");
    }

    // reading real part
    if (func_data.real_size > 0) {
        std::stringstream fname;
        fname << file << "_re";
        alloc(NUMBER::Real, mra);
        real().loadTree(fname.str());
    }

    // reading imaginary part
    if (func_data.imag_size > 0) {
        std::stringstream fname;
        fname << file << "_im";
        alloc(NUMBER::Imag, mra);
        imag().loadTree(fname.str());
    }
    delete mra;
}

/** @brief Returns a character representing the spin (a/b/p) */
char Orbital::printSpin() const {
    char sp = 'u';
    if (this->spin() == SPIN::Paired) sp = 'p';
    if (this->spin() == SPIN::Alpha) sp = 'a';
    if (this->spin() == SPIN::Beta) sp = 'b';
    return sp;
}

} // namespace mrchem
