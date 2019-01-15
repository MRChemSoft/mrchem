/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
#include "PeriodicTable.h"
#include "Element.h"
#include "MRCPP/Printer" //testing only, remove when done

namespace mrchem {

Cavity::Cavity(std::vector<mrcpp::Coord<3>> &coord, std::vector<double> &R, double slope, double eps_i, double eps_o) {
    this->pos = coord;
    this->R = R;
    this->alpha.assign(R.size(), 1.0);
    this->d = slope;
    this->e_i = eps_i;
    this->e_o = eps_o;
}

Cavity::Cavity(const std::vector<std::string> &coord_str,
               double slope,
               double eps_i,
               double eps_o,
               bool atom_based_cavity) {
    this->d = slope;
    this->e_i = eps_i;
    this->e_o = eps_o;
    this->abc = atom_based_cavity;
    readCoordinateString(coord_str);
}

void Cavity::implementEpsilon(bool isinv, bool islinear) {
    this->is_inv = isinv;
    this->is_linear = islinear;
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

    int I = is_inv, L = is_linear;
    I = I << 1;
    int expression = I | L;

    switch (expression) {
        case 0:
            val = e_i * std::exp(log(e_o / e_i) * (1 - C));
            break;
        case 1:
            val = e_i * C + e_o * (1 - C);
            break;
        case 2:
            val = (1 / e_i) * std::exp(log(e_i / e_o) * (1 - C));
            break;
        case 3:
            val = 1 / (e_o + C * (e_i - e_o));
            break;
        default:
            std::cout << "no expression" << std::endl;
    }

    return val;
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

} //namespace mrchem
