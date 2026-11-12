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

#include "utils/gto_utils/DaltonAOReader.h"

#include <MRCPP/Printer>

#include <filesystem>

namespace mrchem {
namespace gto_utils {

DaltonAOReader::DaltonAOReader(const std::string &basisFile) {
    if (not std::filesystem::is_regular_file(basisFile))
        MSG_ABORT("Failed to open basis set file: " << basisFile);
    this->AOFile = basisFile;
    readBasisFile();
}

void DaltonAOReader::readBasisFile() {
    std::ifstream ifs(this->AOFile.c_str());
    size_t nTypes;
    std::string line;
    // First line is a comment
    if (not std::getline(ifs, line))
        MSG_ABORT("Unexpected end of file " << this->AOFile << "  while reading basis sets!");
    // Second line is number of atom types
    if (not std::getline(ifs, line))
        MSG_ABORT("Unexpected end of file " << this->AOFile << "  while reading basis sets!");
    std::istringstream iss(line);
    iss >> nTypes;
    for (size_t i = 0; i < nTypes; i++)
        readAtomBlock(ifs);
}

void DaltonAOReader::readAtomBlock(std::ifstream &ifs) {
    double z;
    size_t nAtoms, nShell;

    ifs >> z;
    ifs >> nAtoms;
    ifs >> nShell;

    IntVector funcsPerShell(nShell);
    for (size_t i = 0; i < nShell; i++)
        ifs >> funcsPerShell(i);

    readAtomData(ifs, nAtoms, z);
    gto_utils::AOBasis bas;
    for (size_t i = 0; i < nShell; i++) {
        for (size_t j = 0; j < funcsPerShell(i); j++)
            readContractionBlock(ifs, bas, i);
    }
    for (size_t i = 0; i < nAtoms; i++)
        this->basis.push_back(std::make_shared<AOBasis>(bas));
}

void DaltonAOReader::readAtomData(std::ifstream &ifs, size_t n_atoms, double z) {
    mrcpp::Coord<3> coord;
    std::string sym;
    for (size_t j = 0; j < n_atoms; j++) {
        ifs >> sym;
        if (sym.size() > 2)
            sym = sym.erase(2);
        for (int d = 0; d < 3; d++)
            ifs >> coord[d];
        PeriodicTable pt;
        const Element &element = pt.getElement(sym.c_str());
        this->nuclei.push_back(std::make_shared<Nucleus>(element, z, coord));
    }
}

void DaltonAOReader::readContractionBlock(std::ifstream &ifs, AOBasis &bas, size_t l) {
    size_t nprim, nctr;
    ifs >> nprim;
    ifs >> nctr;

    size_t start = bas.size();
    for (size_t i = 0; i < nctr; i++) {
        AOContraction ctr(l);
        bas.append(ctr);
    }

    double expo, coef;
    for (size_t k = 0; k < nprim; k++) {
        ifs >> expo;
        for (size_t i = 0; i < nctr; i++) {
            ifs >> coef;
            if (std::abs(coef) > mrcpp::MachineZero)
                bas.getContraction(start + i).append(expo, coef);
        }
    }
}

} // namespace gto_utils
} // namespace mrchem