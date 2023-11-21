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
#include "StepFunction.h"

#include <MRCPP/MWFunctions>

#include "Cavity.h"

namespace mrchem {
namespace detail {
void print_header(const std::string &header, const std::string &formulation, double in_value, double out_value) {
    mrcpp::print::header(0, header);
    print_utils::text(0, "Formulation", formulation, true);
    print_utils::scalar(0, "Value inside Cavity", in_value, "(in)", 6);
    print_utils::scalar(0, "Value outside Cavity", out_value, "(out)", 6);
}
} // namespace detail

StepFunction::StepFunction(std::shared_ptr<mrchem::Cavity> cavity, double val_in, double val_out)
        : in(val_in)
        , out(val_out)
        , cavity{std::move(cavity)} {}

void StepFunction::printParameters() const {
    // Collect relevant quantities
    auto c_pin = this->cavity;

    auto coords = c_pin->getCoordinates();
    auto radii = c_pin->getRadii();
    auto radii_0 = c_pin->getOriginalRadii();
    auto alphas = c_pin->getRadiiScalings();
    auto sigmas = c_pin->getWidths();
    auto betas = c_pin->getWidthScalings();

    // Set widths
    auto w0 = mrcpp::Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 9;
    auto w3 = 6;
    auto w4 = 10;
    auto w5 = w0 - w1 - w2 - 3 * w3 - 3 * w4;

    // Build table column headers
    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "R_0";
    o_head << std::setw(w3 + 1) << "Alpha";
    o_head << std::setw(w3 - 1) << "Beta";
    o_head << std::setw(w3) << "Sigma";
    o_head << std::setw(w5) << "Radius";
    o_head << std::setw(w4) << "x";
    o_head << std::setw(w4) << "y";
    o_head << std::setw(w4) << "z";

    // Print
    printHeader();
    mrcpp::print::separator(0, '-');
    println(0, o_head.str());
    mrcpp::print::separator(0, '-');
    for (auto i = 0; i < coords.size(); i++) {
        auto coord = coords[i];
        auto x = coord[0];
        auto y = coord[1];
        auto z = coord[2];
        auto r = radii[i];
        auto r_0 = radii_0[i];
        auto alpha = alphas[i];
        auto beta = betas[i];
        auto sigma = sigmas[i];

        std::stringstream o_coord;
        o_coord << std::setw(w1) << i;
        o_coord << std::setw(w2) << std::setprecision(4) << std::fixed << r_0;
        o_coord << std::setw(w3) << std::setprecision(2) << std::fixed << alpha;
        o_coord << std::setw(w3) << std::setprecision(2) << std::fixed << beta;
        o_coord << std::setw(w3) << std::setprecision(2) << std::fixed << sigma << "  ->";
        o_coord << std::setw(w5 - 4) << std::setprecision(4) << std::fixed << r;
        o_coord << std::setw(w4) << std::setprecision(6) << std::fixed << x;
        o_coord << std::setw(w4) << std::setprecision(6) << std::fixed << y;
        o_coord << std::setw(w4) << std::setprecision(6) << std::fixed << z;
        println(0, o_coord.str());
    }
    mrcpp::print::separator(0, '=', 2);
}

} // namespace mrchem
