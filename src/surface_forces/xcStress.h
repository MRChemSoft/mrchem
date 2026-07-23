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

#include "mrchem.h"
#include "mrdft/MRDFT.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include <Eigen/Core>
#include <vector>

namespace surface_force {

std::vector<Eigen::Matrix3d> xcLDAStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid);
std::vector<Eigen::Matrix3d> xcLDASpinStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta);
std::vector<Eigen::Matrix3d> xcGGAStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid, Eigen::MatrixXd &nablaRhoGrid);
std::vector<Eigen::Matrix3d> xcGGASpinStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p,
                                             Eigen::MatrixXd &rhoGridAlpha,
                                             Eigen::MatrixXd &rhoGridBeta,
                                             Eigen::MatrixXd &nablaRhoGridAlpha,
                                             Eigen::MatrixXd &nablaRhoGridBeta);
std::vector<Eigen::Matrix3d> getXCStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p,
                                         mrcpp::FunctionTreeVector<3> &xc_pots,
                                         std::shared_ptr<mrchem::OrbitalVector> phi,
                                         std::shared_ptr<mrchem::NablaOperator> nabla,
                                         Eigen::MatrixXd &gridPos,
                                         bool isOpenShell,
                                         double prec);

} // namespace surface_force
