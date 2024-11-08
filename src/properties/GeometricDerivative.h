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

#include <nlohmann/json.hpp>

#include "mrchem.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

class GeometricDerivative final {
public:
    explicit GeometricDerivative(int k = 1)
            : nuclear(math_utils::init_nan(k, 3))
            , electronic(math_utils::init_nan(k, 3)) {}

    DoubleMatrix getTensor() const { return this->nuclear + this->electronic; }
    DoubleMatrix &getNuclear() { return this->nuclear; }
    DoubleMatrix &getElectronic() { return this->electronic; }
    const DoubleMatrix &getNuclear() const { return this->nuclear; }
    const DoubleMatrix &getElectronic() const { return this->electronic; }

    void print(const std::string &id) const {
        mrcpp::print::header(0, "Geometric Derivative (" + id + ")");
        print_utils::matrix(0, "Total", getTensor());
        print_utils::scalar(0, "Norm", getTensor().norm(), "(au)");
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Nuclear", getNuclear());
        print_utils::scalar(0, "Norm", getNuclear().norm(), "(au)");
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Electronic", getElectronic());
        print_utils::scalar(0, "Norm", getElectronic().norm(), "(au)");
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "Est. var. of force error", variance(), "(au)", 4, true);
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        return {{"total", print_utils::eigen_to_vector(getTensor(), 1.0e-12)},
                {"total_norm", getTensor().norm()},
                {"nuclear", print_utils::eigen_to_vector(getNuclear(), 1.0e-12)},
                {"nuclear_norm", getNuclear().norm()},
                {"electronic", print_utils::eigen_to_vector(getElectronic(), 1.0e-12)},
                {"electronic_norm", getElectronic().norm()}};
    }

    /**
     * @brief Compute the variance of the total force using the formula
     * sigma = sqrt(1/(3N) * sum_i sum_j (F_{ij})^2) presented in
     * Gubler et al.: Journal of Computational Physics: X Volume 17, November 2023, 100131
    */
    double variance() const {
        DoubleMatrix forces = getTensor();
        Eigen::Vector3d total_force = forces.colwise().sum();
        double force_variance = total_force.norm() * std::sqrt( 1.0 / (3.0 * forces.rows()) );
        return force_variance;
    }

protected:
    DoubleMatrix nuclear;
    DoubleMatrix electronic;
};

} // namespace mrchem
