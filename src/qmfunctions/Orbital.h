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

#pragma once

#include "QMFunction.h"

/** @class Orbital
 *
 * @brief General complex-valued function with spin
 *
 * Inherits the general features of a complex function from QMFunction which
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

/* POD struct for orbital meta data. Used for simple MPI communication. */
struct OrbitalData {
    int rank_id;
    int spin;
    int occ;
    double error;
};

class Orbital final : public QMFunction {
public:
    Orbital();
    Orbital(int spin, int occ = -1, int rank = -1);

    Orbital(const Orbital &orb);
    Orbital &operator=(const Orbital &orb);
    Orbital paramCopy() const;
    Orbital deepCopy();
    Orbital dagger() const;

    void setOcc(int occ) { this->orb_data.occ = occ; }
    void setSpin(int spin) { this->orb_data.spin = spin; }
    void setRankId(int rank) { this->orb_data.rank_id = rank; }
    void setError(double error) { this->orb_data.error = error; }

    int occ() const { return this->orb_data.occ; }
    int spin() const { return this->orb_data.spin; }
    int rankID() const { return this->orb_data.rank_id; }
    double error() const { return this->orb_data.error; }
    OrbitalData &getOrbitalData() { return this->orb_data; }

    void normalize() { this->rescale(1.0/this->norm()); }
    void orthogonalize(Orbital inp);
    void orthogonalize(OrbitalVector inp_vec);

    void saveOrbital(const std::string &file);
    void loadOrbital(const std::string &file);

    char printSpin() const;
    friend std::ostream& operator<<(std::ostream &o, Orbital orb) { return orb.print(o); }

protected:
    OrbitalData orb_data;

    std::ostream& print(std::ostream &o) const;
};

} //namespace mrchem
