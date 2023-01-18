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

#include "CUBEfunction.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <Eigen/Core>
#include <Eigen/StdVector>

namespace mrchem {

CUBEfunction::CUBEfunction(const int N_atoms,
                           const int N_vals,
                           const std::array<int, 3> N_steps,
                           const mrcpp::Coord<3> origin,
                           const std::array<mrcpp::Coord<3>, 3> Voxel_axes,
                           std::vector<int> Z_n,
                           std::vector<double> cube,
                           std::vector<double> atom_charges,
                           std::vector<mrcpp::Coord<3>> atom_coords)
        : N_atoms(N_atoms)
        , N_val(N_vals)
        , N_steps(N_steps)
        , corner(origin)
        , voxel_axes(Voxel_axes)
        , atom_numbers(Z_n)
        , CUBE(cube)
        , atom_charges(atom_charges)
        , atom_coords(atom_coords) {
    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> voxel_axes(&Voxel_axes[0][0]);
    inv_basis = voxel_axes.transpose().inverse();
}

// Do a quadrature of the file
double CUBEfunction::evalf(const mrcpp::Coord<3> &r) const {

    /** I'm assuming that we can describe any arbitrary point as r = (c_i + d_i)X + (c_j + dc_j)Y + (c_k + dc_k)Z
     * where X, Y, and Z are the three voxel basis vectors descibing the step size and direction of the three main vertexes of the CUBE.
     * Here I am assuming two things, first that the X, Y, and Z basis vectors are orthogonal to each other
     * (X \cdot Y = X \cdot Z=Y \cdot Z = 0),and that they are respectivelly orthogonal to the cartesian basis vectors, for example, with
     * the \hat{x} basis vector:\hat{x} \cdot Y = \hat{x} \cdot Z = 0 and \hat{x} \cdot X = X_x.
     * I believe that this might be generalizable to non-orthogonal basis, but I will not confirm that right now.
     * I also do not want to convert to cartesian coordinates either, as I think it is a waste of time.
     * I then use these d_i, d_j and d_k coefficients as parameters in a trilinear interpolation.
     **/

    double c;
    // perform NX_j \cdot r to find the indices i, j and k of the cubefile.
    Eigen::Map<const Eigen::Vector3d> r_vec(&r[0]);
    Eigen::Map<const Eigen::Vector3d> origin(&corner[0]);
    Eigen::Vector3d coeff = ((inv_basis * (r_vec - origin))); // coefficients i, j and k in r = i*X + j*Y + k*Z assuming basis is orthogonal

    // do the trilinear interpolation naively without loops or any logic (just plug in the equations)
    // Do a sanity check on the point we are evaluating at is in the cube or not. If not return 0.
    // this should be potentially done before any operations, to save time
    std::vector<double> on_edge;
    std::transform(coeff.cbegin(), coeff.cend(), N_steps.cbegin(), std::back_inserter(on_edge), [](const auto &ci, const auto &Ni) {
        double di = ci - (Ni - 1.0);
        return ((di <= -1.0 * (Ni - 1.0)) || (di >= 0.0));
    });

    if (std::any_of(on_edge.cbegin(), on_edge.cend(), [](const auto &ci) { return ci; })) {
        c = 0.0;
    } else {

        auto idx0 = std::floor(coeff(0));
        auto idx1 = std::floor(coeff(1));
        auto idx2 = std::floor(coeff(2));

        auto d_idx0 = coeff(0) - idx0;
        auto d_idx1 = coeff(1) - idx1;
        auto d_idx2 = coeff(2) - idx2;

        auto N_steps1 = N_steps.at(1);
        auto N_steps2 = N_steps.at(2);

        auto c000 = CUBE[(idx0)*N_steps1 * N_steps2 + (idx1)*N_steps2 + (idx2)];
        auto c001 = CUBE[(idx0)*N_steps1 * N_steps2 + (idx1)*N_steps2 + (1 + idx2)];
        auto c010 = CUBE[(idx0)*N_steps1 * N_steps2 + (1 + idx1) * N_steps2 + (idx2)];
        auto c011 = CUBE[(idx0)*N_steps1 * N_steps2 + (1 + idx1) * N_steps2 + (1 + idx2)];
        auto c100 = CUBE[(1 + idx0) * N_steps1 * N_steps2 + (idx1)*N_steps2 + (idx2)];
        auto c101 = CUBE[(1 + idx0) * N_steps1 * N_steps2 + (idx1)*N_steps2 + (1 + idx2)];
        auto c110 = CUBE[(1 + idx0) * N_steps1 * N_steps2 + (1 + idx1) * N_steps2 + (idx2)];
        auto c111 = CUBE[(1 + idx0) * N_steps1 * N_steps2 + (1 + idx1) * N_steps2 + (1 + idx2)];

        auto c00 = c000 * (1 - d_idx0) + c100 * d_idx0;
        auto c01 = c001 * (1 - d_idx0) + c101 * d_idx0;
        auto c10 = c010 * (1 - d_idx0) + c110 * d_idx0;
        auto c11 = c011 * (1 - d_idx0) + c111 * d_idx0;

        auto c0 = c00 * (1 - d_idx1) + c10 * d_idx1;
        auto c1 = c01 * (1 - d_idx1) + c11 * d_idx1;

        c = c0 * (1 - d_idx2) + c1 * d_idx2;
    }

    return c;
}

} // namespace mrchem
