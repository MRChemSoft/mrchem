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

#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include "properties/hirshfeld/HirshfeldInterpolator.h"
#include <Eigen/Dense>
#include <mrchem.h>
#include <string>

/**
 * @brief Class for computing the Hirshfeld partitioning of a molecule. Reads and interpolates the Hirshfeld partitioning data.
 * can create MW representation of the Hirshfeld partitioning functions.
 */
class HirshfeldPartition {

public:
    /**
     * @brief Construct a new Hirshfeld Partition object
     * @param mol The molecule for which the Hirshfeld partitioning is to be computed
     * @param data_dir The directory containing the Hirshfeld partitioning data
     */
    HirshfeldPartition(const mrchem::Molecule &mol, std::string data_dir);

    /**
     * @brief Get the integral rho * w_i for a given atom i
     */
    double getHirshfeldPartitionIntegral(int index, mrcpp::CompFunction<3> &rho, double prec) const;

protected:
    /**
     * @brief Evaluate the analytic, interpolated Hirshfeld partitioning function at a given point
     */
    double evalf(const mrcpp::Coord<3> &r, int iAt) const;

    /**
     * @brief Evaluate the log sum of the exponential of the log atomic densities at a given point.
     */
    double lseLogDens(const mrcpp::Coord<3> &r) const;

    /**
     * @brief The nuclei of the molecule
     */
    std::shared_ptr<mrchem::Nuclei> nucs{nullptr};

    /**
     * @brief The number of nuclei in the molecule
     */
    int nNucs;

    /**
     * @brief The atomic density interpolators for the nuclei
     */
    std::vector<HirshfeldRadInterpolater> logDensities;
};
