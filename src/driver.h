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

#include <nlohmann/json.hpp>

namespace mrchem {

class Molecule;
class CUBEfunction;
namespace driver {

void init_molecule(const nlohmann::json &input, Molecule &mol);
nlohmann::json print_properties(const Molecule &mol);
std::vector<mrchem::CUBEfunction> getCUBEFunction(const nlohmann::json &json_inp);

namespace scf {
nlohmann::json run(const nlohmann::json &input, Molecule &mol);
}
namespace rsp {
nlohmann::json run(const nlohmann::json &input, Molecule &mol);
}

} // namespace driver
} // namespace mrchem
