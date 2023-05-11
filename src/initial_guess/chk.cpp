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

#include <MRCPP/Printer>

#include "chk.h"

#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

bool initial_guess::chk::setup(OrbitalVector &Phi, const std::string &chk_file) {
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation     ", "Compute initial orbitals");
    print_utils::text(0, "Method          ", "Read checkpoint file");
    print_utils::text(0, "Checkpoint file ", chk_file);
    mrcpp::print::separator(0, '~', 2);

    auto success = false;
    auto Psi = orbital::load_orbitals(chk_file);
    if (Psi.size() > 0) {
        success = orbital::compare(Psi, Phi);
        Phi = Psi;
    }
    return success;
}

} // namespace mrchem
