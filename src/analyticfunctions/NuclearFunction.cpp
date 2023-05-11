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

#include "NuclearFunction.h"
#include "chemistry/Nucleus.h"
#include "utils/math_utils.h"
#include "utils/print_utils.h"
#include <MRCPP/Printer>

using mrcpp::Printer;

namespace mrchem {

NuclearFunction::NuclearFunction(const Nuclei &nucs, double smooth_prec, double prec) {
    int pprec = Printer::getPrecision();
    int w0 = Printer::getWidth() - 1;
    int w1 = 5;
    int w2 = 8;
    int w3 = 2 * w0 / 9;
    int w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::string(w4, ' ');
    o_head << std::setw(w3) << "Charge";
    o_head << std::setw(w3) << "Precision";
    o_head << std::setw(w3) << "Smoothing";

    println(2, o_head.str());
    mrcpp::print::separator(2, '-');
    this->prec = prec;

    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double c = detail::nuclear_potential_smoothing(smooth_prec, Z);
        this->push_back(nuc, c);

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << nuc.getElement().getSymbol();
        o_row << std::string(w4, ' ');
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << Z;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth_prec;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << c;
        println(2, o_row.str());
    }
}

void NuclearFunction::push_back(const std::string &atom, const mrcpp::Coord<3> &r, double c) {
    PeriodicTable pt;
    Nucleus nuc(pt.getElement(atom.c_str()), r);
    push_back(nuc, c);
}

void NuclearFunction::push_back(const Nucleus &nuc, double c) {
    this->nuclei.push_back(nuc);
    this->smooth.push_back(c);
    double Z = nuc.getCharge();
    double minPot = -Z * 23 / (c * 3.0 * mrcpp::root_pi); // compatibility with older definition
    this->minPot.push_back(minPot);
}

double NuclearFunction::evalf(const mrcpp::Coord<3> &r) const {
    double result = 0.0;
    // all three method have same lowest value at r = 0
    for (int i = 0; i < this->nuclei.size(); i++) {
        double Z_i = this->nuclei[i].getCharge();
        const mrcpp::Coord<3> &R = this->nuclei[i].getCoord();
        double R1 = math_utils::calc_distance(R, r);
        if (smooth_method == "parabola") {
            // second order, the value and first derivative are equal at R0
            double a = -this->minPot[i];
            double R0 = 1.5 * Z_i / a;
            double b = 0.5 * Z_i / (R0 * R0 * R0);
            if (R1 < R0)
                result += -a + b * R1 * R1;
            else
                result += -Z_i / R1;
        } else if (smooth_method == "minimum") {
            // zero order, just take constant
            result += std::max(-Z_i / R1, this->minPot[i]);
        } else if (smooth_method == "very_smooth") {
            // expensive smoothing that conserve moments
            double c = -1.0 / (3.0 * mrcpp::root_pi);
            double S_i = this->smooth[i];
            R1 /= S_i;
            double partResult = -std::erf(R1) / R1 + c * (std::exp(-R1 * R1) + 16.0 * std::exp(-4.0 * R1 * R1));
            result += Z_i * partResult / S_i;
        } else
            MSG_ABORT("smoothing method not recognized");
    }
    return result;
}

bool NuclearFunction::isVisibleAtScale(int scale, int nQuadPts) const {
    double minSmooth = 1.0;
    if (this->smooth.size() > 0) minSmooth = *std::min_element(this->smooth.begin(), this->smooth.end());
    double stdDeviation = std::pow(minSmooth, -0.5);
    auto visibleScale = static_cast<int>(std::floor(std::log2(nQuadPts * 5.0 * stdDeviation)));
    if (scale < visibleScale) {
        return false;
    } else {
        return true;
    }
}

bool NuclearFunction::isZeroOnInterval(const double *a, const double *b) const {
    int totSplit = 0;
    for (int i = 0; i < this->nuclei.size(); i++) {
        const mrcpp::Coord<3> &R = this->nuclei[i].getCoord();
        int split = 1;
        if (a[0] > R[0] or b[0] < R[0]) split = 0;
        if (a[1] > R[1] or b[1] < R[1]) split = 0;
        if (a[2] > R[2] or b[2] < R[2]) split = 0;
        totSplit += split;
    }
    if (totSplit == 0) {
        return true;
    } else {
        return false;
    }
}

} // namespace mrchem
