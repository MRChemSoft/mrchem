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

#include <string>

#include "qmfunctions/qmfunction_fwd.h"

/** @file gto.h
 *
 * @brief Module for generating initial guess of GTO molecular orbitals
 *
 * The initial_guess::gto namespace provides functionality to setup an
 * initial guess of GTO molecular orbitals. The initial guess requires
 * external input to provide basis set information and MO coefficients.
 */

namespace mrchem {
class Molecule;
class Nucleus;

namespace gto_utils {
class OrbitalExp;
}

namespace initial_guess {
namespace gto {

OrbitalVector setup(double prec,
                    const Molecule &mol,
                    const std::string &bas_file,
                    const std::string &mo_file);
OrbitalVector setup(double prec,
                    const Molecule &mol,
                    const std::string &bas_file,
                    const std::string &moa_file,
                    const std::string &mob_file);
OrbitalVector project_ao(double prec,
                         const std::string &bas_file,
                         int spin,
                         int N = -1);
OrbitalVector project_mo(double prec,
                         const std::string &bas_file,
                         const std::string &mo_file,
                         int spin,
                         int N = -1);
 mrcpp::FunctionTree<3>* project_density(double prec,
                         const Nucleus &nuc,
                         const std::string &bas_file,
                         const std::string &dens_file);

} //namespace gto
} //namespace initial_guess
} //namespace mrchem
