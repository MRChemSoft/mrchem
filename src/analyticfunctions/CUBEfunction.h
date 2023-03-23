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

#include <MRCPP/MWFunctions>

namespace mrchem {

class CUBEfunction final : public mrcpp::RepresentableFunction<3> {
public:
    /**
     * @param[in] N_steps size 3 array of the number of steps in each voxel axis. 0 is the X_axis, 1 is the Y_axis and 2 is the Z_axis
     * @param[in] origin
     * @param[in] Voxel_axes size 3x3 matrix of the voxel axes, first index denotes which voxel, second denotes stepsize on each cartesian coordinate
     * @param[in] cube flattened volumetric data. Indexing here works as  [x_step*N_steps[1]*N_steps[2] + y_step*N_steps[2] + z_step].
     */
    CUBEfunction(const Eigen::Vector3i &N_steps, const Eigen::Vector3d &origin, const Eigen::Matrix3d &Voxel_axes, const Eigen::VectorXd &cube)
            : N_steps(N_steps)
            , origin(origin)
            , inv_basis(Voxel_axes.inverse())
            , CUBE(cube) {}

    /**
     * @param[in] N_steps size 3 array of the number of steps in each voxel axis. 0 is the X_axis, 1 is the Y_axis and 2 is the Z_axis
     * @param[in] origin
     * @param[in] Voxel_axes size 3x3 matrix of the voxel axes, first index denotes which voxel, second denotes stepsize on each cartesian coordinate
     * @param[in] cube flattened volumetric data. Indexing here works as  [x_step*N_steps[1]*N_steps[2] + y_step*N_steps[2] + z_step].
     */
    CUBEfunction(const std::array<int, 3> &N_steps, const mrcpp::Coord<3> &origin, const std::array<mrcpp::Coord<3>, 3> &Voxel_axes, const std::vector<double> &cube)
            : N_steps(N_steps.data())
            , origin(origin.data())
            , CUBE(Eigen::Map<const Eigen::VectorXd>(cube.data(), cube.size())) {
        Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> voxel_axes(&Voxel_axes[0][0]);
        inv_basis = voxel_axes.transpose().inverse();
    }

    double evalf(const mrcpp::Coord<3> &r) const override;

protected:
    Eigen::Vector3i N_steps;
    Eigen::Vector3d origin;
    Eigen::Matrix3d inv_basis;
    Eigen::VectorXd CUBE;
};
} // namespace mrchem
