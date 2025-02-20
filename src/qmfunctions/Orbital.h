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

#include "MRCPP/MWFunctions"
#include "MRCPP/Parallel"
#include "mrchem.h"

/** @class Orbital
 *
 * @brief General complex-valued function with spin
 *
 * Inherits the general features of a complex function from mrcpp::CompFunction<3> which
 * means separate MW function representations for the real and imaginary parts.
 * Note that there are several options for copying/assignment: the proper copy
 * constructor and assignment operator are *shallow* copies, which means that
 * they simply copy the *re and *im pointers (no transfer of ownership).
 * Additionaly, there is a deepCopy() which returns a *full* copy of the orbital,
 * and a paramCopy() which returns an empty orbital with the same rank_id/spin/occ.
 *
 * NOTE: since the standard copies are shallow copies and several orbitals can
 * point to the same MW functions, it is YOUR responibility to keep track of the
 * ownership and free the FunctionTree pointers before the last orbital goes out
 * of scope.
 */

namespace mrchem {

// Note: cannot only define "getSpin()", because sometime we only have a CompFunction, not an Orbital
#define spin() func_ptr->data.n1[0]
#define occ() func_ptr->data.d1[0]
class Orbital : public mrcpp::CompFunction<3> {
public:
    Orbital();
    Orbital(SPIN::type spin);
    Orbital(Orbital &orb);
    Orbital(const mrcpp::CompFunction<3> &orb);
    Orbital(int spin, double occ, int rank = -1);
    Orbital dagger() const;

    //    const int spin() const {return data().n1[0];}
    //    const int occ() const {return data().d1[0];}
    char printSpin() const;
    void setSpin(int spin) {this->func_ptr->data.n1[0] = spin;}
    void saveOrbital(const std::string &file);
    void loadOrbital(const std::string &file);
};

// All MPI processes have a vector of full length, but
// only "my_func" are fully defined.
// The others orbitals (not my_func) have only basic data (spin etc) but no trees defined.
// Vector of orbitals where all orbitals are defined for all MPI, should not have typoe
// OrbitalVectors, but directly vector<Orbital>.
class OrbitalVector : public mrcpp::CompFunctionVector {
public:
    OrbitalVector(int N = 0) : mrcpp::CompFunctionVector(N) {}
    void push_back(Orbital orb) {
        mrcpp::CompFunction<3>& compfunc = orb;
        compfunc.func_ptr->rank = size();
        mrcpp::CompFunctionVector* compfuncvec = this;
        compfuncvec->push_back(compfunc); // we must push in the vector<CompFunction>, not into the OrbitalVector!
    }
    // Overloaded operator[] to return an Orbital element
    // for read (returns lvalue)
    Orbital operator[](int i) const {
        // Create a temporary copy of the base class element to avoid discarding qualifiers
        mrcpp::CompFunction<3> func = mrcpp::CompFunctionVector::operator[](i);
        return Orbital(func);
    }
    // Non-const version of operator[] to allow modification of elements
    // for write (returns rvalue). Cannot return an Orbital
    mrcpp::CompFunction<3>& operator[](int i) {
        return mrcpp::CompFunctionVector::operator[](i);
    }
};

} // namespace mrchem
