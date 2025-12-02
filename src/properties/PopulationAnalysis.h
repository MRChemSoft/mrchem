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
#include "utils/print_utils.h"

#include <nlohmann/json.hpp>
#include <string>

namespace mrchem {
/**
 * @class PopulationAnalysis
 * @brief Store and print population analysis results
 */
class PopulationAnalysis final {
public:
    DoubleMatrix getMatrix() const { return this->populations; } ///< @return Population values

    /**
     * @brief Set population value matrix
     * @param v: population value matrix
     */
    void setMatrix(const DoubleMatrix &v) { this->populations = v; }

    /**
     * @brief Print population analysis
     * @param id: identifier string
     */
    void print(const std::string &id) const {
        auto w0 = mrcpp::Printer::getWidth() - 1;
        auto w1 = 13;
        auto w3 = 2 * w0 / 9;
        auto w4 = w0 - w1 - 3 * w3;

        std::stringstream o_head;
        o_head << std::setw(w1) << "Orbital  ";
        o_head << std::string(w4 - 1, ' ') << ':';
        
        mrcpp::print::header(0, "Orbital Populations (" + id + ")");
        if (populations.cols() == 1) {
            o_head << std::setw(w3*3) << "total";
            println(0, o_head.str());
            mrcpp::print::separator(0, '-');
            for (int i = 0; i < populations.rows(); i++) {
                std::string text = "   " + std::to_string(i);
                print_utils::scalar(0, text, populations(i, 0));
            }
        }
        else {
            o_head << std::setw(w3) << "lower half";
            o_head << std::setw(w3) << "upper half";
            o_head << std::setw(w3) << "total";
            println(0, o_head.str());
            mrcpp::print::separator(0, '-');
            for (int i = 0; i < populations.rows(); i++) {
                std::string text = "   " + std::to_string(i);
                mrcpp::Coord<3> vec = {populations(i, 0), populations(i, 1), populations(i, 2)};
                print_utils::coord(0, text, vec);
            }
        }
        mrcpp::print::separator(0, '=');
    }
    /// @brief Convert the population values to json
    nlohmann::json json() const { return {{"total", print_utils::eigen_to_vector(getMatrix(), 1.0e-12)}}; }

protected:
    DoubleMatrix populations; ///< Population matrix with cols for the different parts of the space
};

} // namespace mrchem