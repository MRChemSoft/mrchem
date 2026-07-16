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

#include "PolyInterpolator.h"
#include <Eigen/Dense>

namespace interpolation_utils {
double polynomialInterpolate5(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x) {
    double xm2 = x_in(0);
    double xm1 = x_in(1);
    double x00 = x_in(2);
    double xp1 = x_in(3);
    double xp2 = x_in(4);

    double ym2 = y_in(0);
    double ym1 = y_in(1);
    double y00 = y_in(2);
    double yp1 = y_in(3);
    double yp2 = y_in(4);
    return ym2 +
           (x - xm2) *
               ((ym1 - ym2) / (xm1 - xm2) +
                (x - xm1) *
                    (((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2) +
                     (x - x00) *
                         ((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1) +
                          ((x - xp1) *
                           (-((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1)) +
                            (-(((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) + ((-y00 + yp1) / (x00 - xp1) + (yp1 - yp2) / (xp1 - xp2)) / (-x00 + xp2)) / (-xm1 + xp2))) /
                              (-xm2 + xp2))));
}

double polynomialInterpolate5_deriv(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x) {
    double xm2 = x_in(0);
    double xm1 = x_in(1);
    double x00 = x_in(2);
    double xp1 = x_in(3);
    double xp2 = x_in(4);

    double ym2 = y_in(0);
    double ym1 = y_in(1);
    double y00 = y_in(2);
    double yp1 = y_in(3);
    double yp2 = y_in(4);

    return (ym1 - ym2) / (xm1 - xm2) +
           (x - xm1) *
               (((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2) +
                (x - x00) * ((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1) +
                             ((x - xp1) *
                              (-((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1)) +
                               (-(((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) + ((-y00 + yp1) / (x00 - xp1) + (yp1 - yp2) / (xp1 - xp2)) / (-x00 + xp2)) / (-xm1 + xp2))) /
                                 (-xm2 + xp2))) +
           (x - xm2) *
               (((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2) +
                (x - x00) * ((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1) +
                             ((x - xp1) *
                              (-((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1)) +
                               (-(((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) + ((-y00 + yp1) / (x00 - xp1) + (yp1 - yp2) / (xp1 - xp2)) / (-x00 + xp2)) / (-xm1 + xp2))) /
                                 (-xm2 + xp2)) +
                (x - xm1) * ((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1) +
                             ((x - x00) *
                              (-((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1)) +
                               (-(((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) + ((-y00 + yp1) / (x00 - xp1) + (yp1 - yp2) / (xp1 - xp2)) / (-x00 + xp2)) / (-xm1 + xp2))) /
                                 (-xm2 + xp2) +
                             ((x - xp1) *
                              (-((-(((y00 - ym1) / (x00 - xm1) + (-ym1 + ym2) / (xm1 - xm2)) / (x00 - xm2)) + ((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) / (-xm2 + xp1)) +
                               (-(((-y00 + ym1) / (x00 - xm1) + (y00 - yp1) / (x00 - xp1)) / (-xm1 + xp1)) + ((-y00 + yp1) / (x00 - xp1) + (yp1 - yp2) / (xp1 - xp2)) / (-x00 + xp2)) / (-xm1 + xp2))) /
                                 (-xm2 + xp2)));
}

int binarySearch(const Eigen::VectorXd &x, const double &x0) {
    int i = 0;
    int j = x.size() - 1;
    int k;
    while (j - i > 1) {
        k = (i + j) / 2;
        if (x0 < x(k)) {
            j = k;
        } else {
            i = k;
        }
    }
    return i;
}
} // namespace interpolation_utils