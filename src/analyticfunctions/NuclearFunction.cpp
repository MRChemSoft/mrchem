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

/*
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

    for (size_t k = 0; k < nucs.size(); k++) {
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
*/
void NuclearFunction::push_back(const std::string &atom, const mrcpp::Coord<3> &r, double p1, double p2) {
    PeriodicTable pt;
    Nucleus nuc(pt.getElement(atom.c_str()), r);
    push_back(nuc, p1, p2);
}

void NuclearFunction::push_back(const Nucleus &nuc, double p1, double p2) {
    this->nuclei.push_back(nuc);
    this->param1.push_back(p1);
    this->param2.push_back(p2);
}

bool NuclearFunction::isVisibleAtScale(int scale, int nQuadPts) const {
    double stdDeviation = 1.0;
    auto visibleScale = static_cast<int>(std::floor(std::log2(nQuadPts * 5.0 * stdDeviation)));
    return (scale >= visibleScale);
}

bool NuclearFunction::isZeroOnInterval(const double *a, const double *b) const {
    int totSplit = 0;
    for (size_t i = 0; i < this->nuclei.size(); i++) {
        const mrcpp::Coord<3> &R = this->nuclei[i].getCoord();
        int split = 1;
        if (a[0] > R[0] or b[0] < R[0]) split = 0;
        if (a[1] > R[1] or b[1] < R[1]) split = 0;
        if (a[2] > R[2] or b[2] < R[2]) split = 0;
        totSplit += split;
    }
    return (totSplit == 0);
}

} // namespace mrchem
