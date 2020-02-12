/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/Printer"

#include "mw.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

OrbitalVector initial_guess::mw::setup(const std::string &file) {
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Read orbitals from file (MW)");
    print_utils::text(0, "File name   ", file);
    mrcpp::print::separator(0, '~', 2);

    return orbital::load_orbitals(file);
}

} // namespace mrchem
