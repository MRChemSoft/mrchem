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

namespace mrchem {

/** @brief Default constructor
 *
 * Initializes with NULL tree pointers.
 */
Orbital::Orbital()
        : mrcpp::CompFunction<3>() {}

/** @brief Constructor with only spin
 *
 * Initializes with NULL tree pointers.
 */
Orbital::Orbital(SPIN::type spin)
        : mrcpp::CompFunction<3>(spin) {
    if (this->spin() < 0) INVALID_ARG_ABORT;
    // d1 is used to store occupancy
    if (this->spin() == SPIN::Paired) this->func_ptr->data.d1[0] = 2;
    if (this->spin() == SPIN::Alpha) this->func_ptr->data.d1[0] = 1;
    if (this->spin() == SPIN::Beta) this->func_ptr->data.d1[0] = 1;
}

/** @brief Constructor
 *
 * @param spin: electron spin (SPIN::Alpha/Beta/Paired)
 * @param occ: occupation
 * @param rank: position in vector if part of a vector
 *
 */
Orbital::Orbital(int spin, double occ, int rank)
        : mrcpp::CompFunction<3>(spin){
    if (this->spin() < 0) INVALID_ARG_ABORT;
    if (this->occ() < 0) {
        // d1 is defined as occupancy
        if (this->spin() == SPIN::Paired) this->func_ptr->data.d1[0] = 2;
        if (this->spin() == SPIN::Alpha) this->func_ptr->data.d1[0] = 1;
        if (this->spin() == SPIN::Beta) this->func_ptr->data.d1[0] = 1;
    }
    this->func_ptr->rank = rank;
}


/** @brief Copy constructor
 *
 * @param orb: orbital to copy
 *
 * Shallow copy: meta data is copied along with the pointers,
 * NO transfer of ownership:
 * both orbitals are pointing to the same tree
 */
Orbital::Orbital(Orbital &orb)
        : mrcpp::CompFunction<3>(orb) {}


/** @brief Copy constructor
 *
 * @param orb: orbital to copy
 *
 * Shallow copy: meta data is copied along with the pointers,
 * NO transfer of ownership:
 * both orbitals are pointing to the same tree
 */
//Orbital::Orbital(const mrcpp::CompFunction<3> &orb)
//        : mrcpp::CompFunction<3>(orb) {}
Orbital::Orbital(const mrcpp::CompFunction<3>& orb)
    : mrcpp::CompFunction<3>(orb) {}


/** @brief Complex conjugation
 *
 * Returns a new orbital which is a shallow copy of *this orbital, with a flipped
 * conjugate parameter. Pointer ownership is not transferred, so *this and the output
 * orbital points to the same MW representations of the function tree,
 * however, they interpret the imaginary part with opposite sign.
 */
Orbital Orbital::dagger() const {
    Orbital out(*this); // Shallow copy
    out.func_ptr->conj = not this->conjugate();
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
    if (isreal()) CompD[0]->saveTree(file);
    if (iscomplex()) CompC[0]->saveTree(file);
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
    if (isreal()) CompD[0]->loadTree(file);
    if (iscomplex()) CompC[0]->loadTree(file);
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
