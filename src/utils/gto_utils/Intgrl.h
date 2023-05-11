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

/** Class to parse and process basis sets in MOLECULE or INTGRL format.
 */
#pragma once

#include <MRCPP/Gaussians>

#include <fstream>
#include <string>
#include <vector>

namespace mrchem {
class Nucleus;

namespace gto_utils {
class AOBasis;

class Intgrl final {
public:
    Intgrl(const std::string &file);
    ~Intgrl();

    int getNNuclei() const { return this->nuclei.size(); }

    const Nucleus &getNucleus(int i) const { return *this->nuclei[i]; }
    const AOBasis &getAOBasis(int i) const { return *this->basis[i]; }

    Nucleus &getNucleus(int i) { return *this->nuclei[i]; }
    AOBasis &getAOBasis(int i) { return *this->basis[i]; }

    mrcpp::GaussExp<3> getAtomBasis(int i, bool norm = true) const;
    mrcpp::GaussExp<3> getMolBasis(bool norm = true) const;

protected:
    std::vector<Nucleus *> nuclei;
    std::vector<AOBasis *> basis;

    void readIntgrlFile(const std::string &fname, std::ifstream &ifs);
    void readContractionBlock(std::ifstream &ifs, AOBasis &basis, int l);
    void readAtomBlock(std::ifstream &ifs);
    void readAtomData(std::ifstream &ifs, int n_atoms, double z);
};

} // namespace gto_utils
} // namespace mrchem
