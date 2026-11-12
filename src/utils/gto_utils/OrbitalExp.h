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

#include <vector>
#include <memory>

#include "MRCPP/Gaussians"

#include "mrchem.h"

namespace mrchem {

class AOFileReader;
namespace gto_utils {

struct CartToSphTransformation {
    std::vector<std::vector<int>> inds;
    std::vector<std::vector<double>> coeffs;
};

class OrbitalExp final {
public:
    OrbitalExp(AOFileReader &reader);
    ~OrbitalExp() = default;

    /**
     * @brief Get the number of orbitals in the expansion.
     * @return The number of orbitals.
     */
    size_t size() const { return this->orbitals.size(); }
    
    /**
     * @brief Get the i-th atomic orbital in the expansion.
     * @param i The index of the atomic orbital.
     * @return A GaussExp<3> object representing the i-th atomic orbital.
     */
    mrcpp::GaussExp<3> getAO(int i) const { return this->orbitals[i]; }

    /**
     * @brief Get the molecular orbital corresponding to the i-th row of the coefficient matrix M.
     * @param i The index of the molecular orbital.
     * @param M The coefficient matrix.
     * @param threshold The threshold for including contributions from atomic orbitals.
     * @return A GaussExp<3> object representing the i-th molecular orbital.
     */
    mrcpp::GaussExp<3> getMO(int i, const DoubleMatrix &M, const double threshold=mrcpp::MachineZero) const;

    /**
     * @brief Get the density corresponding to the density matrix D.
     * @param D The density matrix.
     * @return A GaussExp<3> object representing the density.
     */
    mrcpp::GaussExp<3> getDens(const DoubleMatrix &D) const;
    
    /**
     * @brief Rotate the orbitals using the unitary matrix U.
     * @param U The unitary matrix for rotation.
     */
    void rotate(const DoubleMatrix &U);
    
private:
    bool cartesian;
    std::vector<mrcpp::GaussExp<3>> orbitals;
    
    std::vector<CartToSphTransformation> sph_transformation_data;
    
    /**
     * @brief Get the angular momentum quantum number of the n-th orbital.
     * @param n The index of the orbital.
     * @return The angular momentum quantum number.
     */
    int getAngularMomentum(int n) const;

    /**
     * @brief Read the atomic orbital expansion from the AOFileReader.
     * @param reader The AOFileReader object to read from.
     */
    void readAOExpansion(AOFileReader &reader);

    /**
     * @brief Transform the orbitals from Cartesian to spherical harmonics representation.
     */
    void transformToSpherical();

    /**
     * @brief Get the transformation data for converting Cartesian to spherical harmonics for a given angular momentum.
     * @param l The angular momentum quantum number.
     * @return A reference to the CartToSphTransformation object containing the transformation data.
     */
    CartToSphTransformation &getSphTransformation(int l);
};

} // namespace gto_utils
} // namespace mrchem
