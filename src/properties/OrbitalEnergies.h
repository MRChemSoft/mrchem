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

/** @class OrbitalEnergies
 *
 * @brief Simple POD container to hold the orbital eigenvalues
 *
 */

namespace mrchem {

// clang-format off
class OrbitalEnergies final {
public:
    IntVector &getSpin() { return this->spin; }
    IntVector &getOccupation() { return this->occupation; }
    DoubleVector &getEpsilon() { return this->epsilon; }

    const IntVector &getSpin() const { return this->spin; }
    const IntVector &getOccupation() const { return this->occupation; }
    const DoubleVector &getEpsilon() const { return this->epsilon; }

    void print(const std::string &id) const {
        auto pprec = 2 * mrcpp::Printer::getPrecision();
        auto w0 = mrcpp::Printer::getWidth() - 1;
        auto w1 = 5;
        auto w2 = 2 * w0 / 9;
        auto w3 = w0 - 3 * w1 - 3 * w2;

        std::stringstream o_head;
        o_head << std::setw(w1) << "n";
        o_head << std::setw(w1) << "Occ";
        o_head << std::setw(w1) << "Spin";
        o_head << std::string(w3 - 1, ' ') << ':';
        o_head << std::setw(3 * w2) << "Epsilon";

        mrcpp::print::header(0, "Orbital Energies (" + id + ")");
        println(0, o_head.str());
        mrcpp::print::separator(0, '-');

        for (int i = 0; i < this->epsilon.size(); i++) {
            auto sp = 'u';
            if (this->spin(i) == SPIN::Paired) sp = 'p';
            if (this->spin(i) == SPIN::Alpha) sp = 'a';
            if (this->spin(i) == SPIN::Beta) sp = 'b';
            std::stringstream o_txt;
            o_txt << std::setw(w1 - 1) << i;
            o_txt << std::setw(w1) << this->occupation(i);
            o_txt << std::setw(w1) << sp;
            print_utils::scalar(0, o_txt.str(), this->epsilon(i), "(au)", pprec);
        }
        auto sum_occupied = this->occupation.cast<double>().dot(this->epsilon);
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "Sum occupied", sum_occupied, "(au)", pprec);
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        DoubleVector occ = getOccupation().cast<double>();
        const DoubleVector &eps = getEpsilon();
        std::vector<std::string> spn;
        for (auto i = 0; i < spin.size(); i++) {
            if (this->spin(i) == SPIN::Paired) spn.push_back("p");
            if (this->spin(i) == SPIN::Alpha) spn.push_back("a");
            if (this->spin(i) == SPIN::Beta) spn.push_back("b");
        }
        return {
            {"spin", spn},
            {"occupation", print_utils::eigen_to_vector(occ, 1.0e-12)},
            {"energy", print_utils::eigen_to_vector(eps, 1.0e-12)},
            {"sum_occupied", occ.dot(eps)}
        };
    }

private:
    IntVector spin;
    IntVector occupation;
    DoubleVector epsilon;
};
// clang-format on

} // namespace mrchem
