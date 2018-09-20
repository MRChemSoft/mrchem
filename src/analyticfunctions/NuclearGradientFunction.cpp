/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "MRCPP/Printer"

#include "NuclearGradientFunction.h"
#include "utils/math_utils.h"

namespace mrchem {

double NuclearGradientFunction::evalf(const double *r) const {
    double c_A = this->smooth;
    double Z_A = this->nucleus.getCharge();
    const double *R_A = this->nucleus.getCoord();
    double x_A =  R_A[this->dir] - r[this->dir];
    double r_A = math_utils::calc_distance(R_A, r);
    return -Z_A*(x_A/r_A)*du_dr(r_A/c_A)/(c_A*c_A);
}

bool NuclearGradientFunction::isVisibleAtScale(int scale, int nQuadPts) const {
    double stdDeviation = std::pow(this->smooth, -0.5);
    int visibleScale = int (std::floor(std::log2(nQuadPts*5.0*stdDeviation)));
    return (scale < visibleScale) ? false : true;
}

bool NuclearGradientFunction::isZeroOnInterval(const double *a, const double *b) const {
    bool out = false;
    const double *R = this->nucleus.getCoord();
    if (a[0] > R[0] or b[0] < R[0]) out = true;
    if (a[1] > R[1] or b[1] < R[1]) out = true;
    if (a[2] > R[2] or b[2] < R[2]) out = true;
    return out;
}

double NuclearGradientFunction::du_dr(double r1) const {
    double out = 0.0;
    if (r1 > 6.0) {
        double r2 = r1*r1;
        out = -1.0/r2;
    } else if (r1 > 0.1) {
        double r2 = r1*r1;
        out =   2.0*std::exp(-r2)/(mrcpp::root_pi*r1)
              - std::erf(r1)/r2
              - 1.0/(3.0*mrcpp::root_pi)
              * (2.0*r1*std::exp(-r2) + 128.0*r1*std::exp(-4.0*r2));
    } else if (r1 > 0.0) {
        double r2 = r1*r1;
        double r3 = r1*r2;
        double r5 = r3*r2;
        double r7 = r5*r2;
        double r9 = r7*r2;
        out = - 4.0/3.0*r1
              + 4.0/5.0*r3
              - 2.0/7.0*r5
              + 2.0/27.0*r7
              - 1.0/66.0*r9
              - 1.0/(3.0*mrcpp::root_pi)
              * (2.0*r1*std::exp(-r2) + 128.0*r1*std::exp(-4.0*r2));
    } else {
        MSG_FATAL("Negative radius");
    }
    return out;
}

} //namespace mrchem
