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
#include "utils/print_utils.h"

/** @class ExcitationEnergies
 *
 * @brief Simple POD container to hold the Excitation energies
 *
 */

namespace mrchem {

// clang-format off
class ExcitationEnergies final {
public:
    std::vector<double> &getOmega() { return this->omega; }

    const std::vector<double> &getOmega() const { return this->omega; }

    void print(const std::string &id) const {
        auto pprec = 2 * mrcpp::Printer::getPrecision();
        auto w0 = mrcpp::Printer::getWidth() - 1;
        auto w1 = 5;
        auto w2 = 2 * w0 / 9;
        auto w3 = w0 - 3 * w1 - 3 * w2;

        std::stringstream o_head;
        o_head << std::setw(w1) << "r";
        o_head << std::string(w1 + w1 + w3 - 1, ' ') << ':';
        o_head << std::setw(3 * w2) << "Omega";

        // could use this to print both the starting guess and final 
        // converged frequencies, could help to check the quality of the guess
        mrcpp::print::header(0, "Excitation Energies (" + id + ")"); 
        println(0, o_head.str());
        mrcpp::print::separator(0, '-');

        for (int i = 0; i < this->omega.size(); i++) {
            std::stringstream o_txt;
            o_txt << std::setw(w1 - 1) << i+1;
            print_utils::scalar(0, o_txt.str(), this->omega[i], "(au)", pprec);
        }
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        const std::vector<double> &ome = getOmega();
        nlohmann::json json_out;
        for (auto i = 0; i < ome.size(); i++) {
            std::stringstream state_key;
            state_key << "state-" << i+1;
            json_out[state_key.str()] = ome[i];
        }
        return json_out;
    }

private:
    std::vector<double> omega;
};
// clang-format on

} // namespace mrchem
