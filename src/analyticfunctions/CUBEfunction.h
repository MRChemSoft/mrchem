/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
#include <MRCPP/MWFunctions>

namespace mrchem {

class CUBEfunction final : public mrcpp::RepresentableFunction<3> {
public:
    CUBEfunction(const int N_atoms,
                 const int N_vals,
                 const std::array<int, 3> N_steps,
                 const mrcpp::Coord<3> origin,
                 const std::array<mrcpp::Coord<3>, 3> Voxel_axes,
                 std::vector<int> Z_n,
                 std::vector<double> cube,
                 std::vector<double> atom_charges,
                 std::vector<mrcpp::Coord<3>> atom_coords);
    double evalf(const mrcpp::Coord<3> &r) const override;

protected:
    void normalize_basis();

    int N_atoms;
    int N_val;
    std::array<int, 3> N_steps; // size 3 array of the number of steps in each voxel axis. 0 is the X_axis, 1 is the Y_axis and 2 is the Z_axis
    mrcpp::Coord<3> corner;
    std::array<mrcpp::Coord<3>, 3> voxel_axes; // size 3x3 matrix of the voxel axes, first index denotes which voxel, second denotes stepsize on each cartesian coordinate
    std::vector<int> atom_numbers;
    std::vector<double> CUBE; // indexing here works as  [x_step*N_steps[1]*N_steps[2] + y_step*N_steps[2] + z_step].
    std::vector<double> atom_charges;
    std::vector<mrcpp::Coord<3>> atom_coords;

    Eigen::Matrix3d normalized_basis; // multiply each row by its 1/norm^2
};
} // namespace mrchem
