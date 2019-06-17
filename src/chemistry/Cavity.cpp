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

#include <array>
#include <cmath>
#include <vector>

#include "Cavity.h"
#include "Element.h"
#include "MRCPP/Printer" //testing only, remove when done
#include "PeriodicTable.h"
#include "utils/math_utils.h"

namespace mrchem {
// Add different dimensions to the cavity

Cavity::Cavity(std::vector<mrcpp::Coord<3>> &coord, std::vector<double> &R, double slope) {
    this->pos = coord;
    this->R = R;
    this->alpha.assign(R.size(), 1.0);
    this->d = slope;
}

Cavity::Cavity(const std::vector<std::string> &coord_str, double slope, bool atom_based_cavity) {
    this->d = slope;
    this->abc = atom_based_cavity;
    readCoordinateString(coord_str);
}

double Cavity::evalf(const mrcpp::Coord<3> &r) const {
    double C = 1.0;
    double s, O;
    double val;
    for (int i = 0; i < pos.size(); i++) {
        s = math_utils::calc_distance(pos[i], r) - R[i] * alpha[i]; // scale vdwr by 1.2
        O = 0.5 * (1 + std::erf(s / d));
        C *= 1 - (1 - O);
    }
    C = 1 - C;
    return C;
}

void Cavity::readCoordinateString(const std::vector<std::string> &coord_str) {
    int nAtoms = coord_str.size();
    mrcpp::Coord<3> coord;
    double Rad;
    std::string sym;
    if (this->abc) {
        PeriodicTable P;
        for (int i = 0; i < nAtoms; i++) {
            std::stringstream ss;
            ss.str(coord_str[i]);
            ss >> sym;
            ss >> coord[0];
            ss >> coord[1];
            ss >> coord[2];
            Rad = P.getElement(sym.c_str()).getVdw();
            this->pos.push_back(coord);
            this->R.push_back(Rad);
        }
        this->alpha.assign(R.size(), 1.2);

    } else if (not this->abc) {
        for (int i = 0; i < nAtoms; i++) {
            std::stringstream ss;
            ss.str(coord_str[i]);
            ss >> sym;
            ss >> coord[0];
            ss >> coord[1];
            ss >> coord[2];
            ss >> Rad;
            this->pos.push_back(coord);
            this->R.push_back(Rad);
        }
        this->alpha.assign(R.size(), 1.0);
    }
}

bool Cavity::isVisibleAtScale(int scale, int nQuadPts) const {

    double stdDeviation = d;
    auto visibleScale = static_cast<int>(-std::floor(std::log2(nQuadPts * 2.0 * stdDeviation)));

    if (scale < visibleScale) { return false; }

    return true;
}
// maybe necessary to state value outside the boundaries of the cavity, that is
// eps_o and eps_i
bool Cavity::isZeroOnInterval(const double *a, const double *b) const {
    for (int k = 0; k < pos.size(); k++) {
        for (int i = 0; i < 3; i++) {
            double stdDeviation = d;
            double cavityMinOut = (this->pos[k][i] - R[i]) - 3.0 * stdDeviation;
            double cavityMinIn = (this->pos[k][i] - R[i]) + 3.0 * stdDeviation;
            double cavityMaxIn = (this->pos[k][i] + R[i]) - 3.0 * stdDeviation;
            double cavityMaxOut = (this->pos[k][i] + R[i]) + 3.0 * stdDeviation;
            if (a[i] > cavityMaxOut or (a[i] > cavityMinIn and b[i] < cavityMaxIn) or b[i] < cavityMinOut) {
                return true;
            }
        }
    }
    return false;
}

} // namespace mrchem
