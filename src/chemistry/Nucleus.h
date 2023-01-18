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

#include <iostream>
#include <vector>

#include "MRCPP/MWFunctions"

#include "Element.h"
#include "PeriodicTable.h"

namespace mrchem {

class Nucleus final {
public:
    Nucleus(const Element &elm, const mrcpp::Coord<3> &r)
            : charge(elm.getZ())
            , coord(r)
            , element(&elm) {}
    Nucleus(const Nucleus &nuc)
            : charge(nuc.charge)
            , coord(nuc.coord)
            , element(nuc.element) {}
    Nucleus &operator=(const Nucleus &nuc) {
        if (this != &nuc) {
            this->charge = nuc.charge;
            this->element = nuc.element;
            this->coord = nuc.coord;
        }
        return *this;
    }

    void setCharge(double z) { this->charge = z; }
    void setCoord(const mrcpp::Coord<3> &r) { this->coord = r; }

    double getCharge() const { return this->charge; }
    const mrcpp::Coord<3> &getCoord() const { return this->coord; }
    const Element &getElement() const { return *this->element; }

    friend std::ostream &operator<<(std::ostream &o, const Nucleus &nuc) {
        o << std::endl << *nuc.element << "   Z = ";
        o << nuc.charge << std::endl << "   ";
        for (int d = 0; d < 3; d++) o << nuc.coord[d] << " ";
        o << std::endl;
        return o;
    }

private:
    double charge;
    mrcpp::Coord<3> coord;
    const Element *element;
};

class Nuclei : public std::vector<Nucleus> {
public:
    void push_back(const Nucleus &nuc) { std::vector<Nucleus>::push_back(nuc); }
    void push_back(const std::string &atom, const mrcpp::Coord<3> &xyz) {
        PeriodicTable pt;
        Nucleus nuc(pt.getElement(atom.c_str()), xyz);
        std::vector<Nucleus>::push_back(nuc);
    }
};

} // namespace mrchem
