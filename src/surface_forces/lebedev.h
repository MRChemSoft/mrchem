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

#include <Eigen/Dense>
#include <string>

class LebedevIntegrator {
public:
    Eigen::MatrixXd points;  /**< Matrix storing the Cartesian coordinates of the Lebedev grid points. */
    Eigen::VectorXd weights; /**< Vector storing the weights associated with each Lebedev grid point. */
    Eigen::MatrixXd normals; /**< Matrix storing the normal vectors (to the integration sphere) at each Lebedev grid point. */
    int n;                   /**< Number of Lebedev grid points. */

    LebedevIntegrator(int nPoints, double radius, const Eigen::Vector3d &center);
    Eigen::MatrixXd getPoints() const;
    Eigen::VectorXd getWeights() const;
    Eigen::MatrixXd getNormals() const;

private:
    void getLebedevData(int nPoints);
    void readLebedevFile(const std::string &filename);
    void calculateCartesianPoints(double radius, const Eigen::Vector3d &center);
};