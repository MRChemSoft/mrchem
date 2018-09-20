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

/** 
 *
 * \date Jun 7, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 *
 */

#pragma once

#include <string>
#include <map>

#include "chemistry_fwd.h"

namespace mrchem {

class PeriodicTable final {
public:
    PeriodicTable() { }
    ~PeriodicTable() { }
    const Element &getElement(int Z) const;
    const Element &getElement(const char *id) const;
    static const int nElements=112;
protected:
    typedef std::map<std::string, const Element *> map_t;
    static const Element elements[];
    static map_t byName; 
    static map_t bySymbol; 
private:
    static map_t _init_byname();
    static map_t _init_bysymbol();
};

} //namespace mrchem
