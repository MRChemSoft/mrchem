/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

// clang-format off
class NMRShielding final {
public:
    explicit NMRShielding(const mrcpp::Coord<3> &k = {}, const mrcpp::Coord<3> &o = {}) : r_K(k), r_O(o) {}

    void setCoordK(const mrcpp::Coord<3> &k) { this->r_K = k; }
    void setOrigin(const mrcpp::Coord<3> &o) { this->r_O = o; }

    const mrcpp::Coord<3> &getCoordK() const { return this->r_K; }
    const mrcpp::Coord<3> &getOrigin() const { return this->r_O; }

    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    void print(const std::string &id) const {
        auto sigma = getTensor();
        DoubleVector diag = math_utils::init_nan(3);
        if (not sigma.hasNaN()) {
            Eigen::EigenSolver<DoubleMatrix> es;
            es.compute(sigma);
            if (es.eigenvalues().imag().norm() > 1.0e-6) MSG_WARN("Complex NMR eigenvalue")

            diag = es.eigenvalues().real();
        }
        auto iso_ppm = diag.sum() / 3.0;
        auto ani_ppm = (3.0/2.0)*diag.maxCoeff() - diag.sum() / 2.0;

        mrcpp::print::header(0, "NMR shielding (" + id + ")");
        print_utils::coord(0, "r_O", getOrigin());
        print_utils::coord(0, "r_K", getCoordK());
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Diamagnetic", getDiamagnetic());
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Paramagnetic", getParamagnetic());
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Total tensor", getTensor());
        mrcpp::print::separator(0, '-');
        print_utils::vector(0, "Diagonalized tensor", diag);
        print_utils::scalar(0, "Isotropic average", iso_ppm, "(ppm)");
        print_utils::scalar(0, "Anisotropy", ani_ppm, "(ppm)");
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        auto sigma = getTensor();
        DoubleVector diag = math_utils::init_nan(3);
        if (not sigma.hasNaN()) {
            Eigen::EigenSolver<DoubleMatrix> es;
            es.compute(sigma);
            diag = es.eigenvalues().real();
        }
        auto iso_ppm = diag.sum() / 3.0;
        auto ani_ppm = (3.0/2.0)*diag.maxCoeff() - diag.sum() / 2.0;

        return {
            {"r_K", getCoordK()},
            {"r_O", getOrigin()},
            {"tensor_dia", print_utils::eigen_to_vector(getDiamagnetic(), 1.0e-12)},
            {"tensor_para", print_utils::eigen_to_vector(getParamagnetic(), 1.0e-12)},
            {"tensor", print_utils::eigen_to_vector(getTensor(), 1.0e-12)},
            {"diagonalized_tensor", print_utils::eigen_to_vector(diag, 1.0e-12)},
            {"isotropic_average", iso_ppm},
            {"anisotropy", ani_ppm}
        };
    }

private:
    mrcpp::Coord<3> r_K;
    mrcpp::Coord<3> r_O;
    DoubleMatrix dia_tensor{math_utils::init_nan(3,3)};
    DoubleMatrix para_tensor{math_utils::init_nan(3,3)};
};
// clang-format on

} // namespace mrchem
