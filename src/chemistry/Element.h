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

#include <string>

/** Basic chemical element data. */
namespace mrchem {

class Element {
public:
    Element(int z, const char *s, const char *n, double m, double vdw, double cov, double g)
            : Z(z)
            , symbol(s)
            , name(n)
            , mass(m)
            , r_vdw(vdw)
            , r_cov(cov)
            , g_val(g) {}
    virtual ~Element() {}

    const std::string &getSymbol() const { return this->symbol; }
    const std::string &getName() const { return this->name; }

    int getZ() const { return this->Z; }
    double getMass() const { return this->mass; }
    double getVdw() const { return this->r_vdw; }
    double getCov() const { return this->r_cov; }
    double getGValue() const {
        if (this->g_val < 0.0) std::abort();
        return this->g_val;
    }

    friend std::ostream &operator<<(std::ostream &o, const Element &e) {
        o << e.symbol;
        return o;
    }

protected:
    const int Z;              /** atomic number */
    const std::string symbol; /** atomic symbol */
    const std::string name;   /** element name */
    const double mass;        /** atomic mass */
    const double r_vdw;       /** van der waals radius */
    const double r_cov;       /** covalent radius */
    const double g_val;       /** g-value */
};

} // namespace mrchem
