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

#include <iomanip>
#include <iostream>

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "print_utils.h"

#include "qmfunctions/QMFunction.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

std::vector<double> print_utils::eigen_to_vector(const DoubleVector &inp, double thrs) {
    std::vector<double> out;
    for (int i = 0; i < inp.size(); i++) {
        auto out_i = inp(i);
        if (std::abs(out_i) < thrs) out_i = 0.0;
        out.push_back(out_i);
    }
    return out;
}

std::vector<double> print_utils::eigen_to_vector(const DoubleMatrix &inp, double thrs) {
    std::vector<double> out;
    for (int i = 0; i < inp.rows(); i++) {
        for (int j = 0; j < inp.cols(); j++) {
            auto out_ij = inp(i, j);
            if (std::abs(out_ij) < thrs) out_ij = 0.0;
            out.push_back(out_ij);
        }
    }
    return out;
}

std::string print_utils::dbl_to_str(double d, int p, bool sci) {
    std::stringstream o_dbl;
    if (sci) {
        o_dbl << std::setprecision(p) << std::scientific << d;
    } else {
        o_dbl << std::setprecision(p) << std::fixed << d;
    }
    return o_dbl.str();
}

void print_utils::headline(int level, const std::string &txt) {
    auto pwidth = Printer::getWidth();
    auto txt_width = txt.size();
    auto pre_spaces = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');

    mrcpp::print::separator(level, ' ');
    mrcpp::print::separator(level, '*');
    println(level, pre_str << std::string(txt_width, ' ') << post_str);
    println(level, pre_str << txt << post_str);
    println(level, pre_str << std::string(txt_width, ' ') << post_str);
    mrcpp::print::separator(level, '*');
    mrcpp::print::separator(level, ' ');
    mrcpp::print::separator(level, ' ');
}

void print_utils::text(int level, const std::string &txt, const std::string &val) {
    int w0 = Printer::getWidth() - 2;
    int w1 = w0 * 2 / 9;
    int w2 = w0 - 3 * w1;
    int w3 = w2 - (txt.size() + 1);

    std::stringstream o;
    o << " " << txt << std::string(w3, ' ') << ": " << val;
    println(level, o.str());
}

void print_utils::coord(int level, const std::string &txt, const mrcpp::Coord<3> &val, int p, bool s) {
    if (p < 0) p = Printer::getPrecision();
    int w0 = Printer::getWidth() - 2;
    int w1 = w0 * 2 / 9;
    int w2 = w0 - 3 * w1;
    int w3 = w2 - (txt.size() + 1);

    std::stringstream o;
    o << " " << txt << std::string(w3, ' ') << ":";
    if (s) {
        for (auto r : val) o << std::setw(w1) << std::setprecision(p) << std::scientific << r;
    } else {
        for (auto r : val) o << std::setw(w1) << std::setprecision(p) << std::fixed << r;
    }
    println(level, o.str());
}

void print_utils::scalar(int level, const std::string &txt, double val, const std::string &unit, int p, bool s) {
    if (p < 0) p = Printer::getPrecision();
    int w0 = Printer::getWidth() - 2;
    int w1 = w0 * 2 / 9;
    int w2 = w0 - 3 * w1;
    int w3 = w2 - (txt.size() + 1);

    std::stringstream o;
    o << " " << txt << std::string(w3, ' ') << ":";
    o << std::setw(w1) << unit;
    if (std::isnan(val)) {
        o << std::setw(2 * w1) << "N/A";
    } else if (s) {
        o << std::setw(2 * w1) << std::setprecision(p) << std::scientific << val;
    } else {
        o << std::setw(2 * w1) << std::setprecision(p) << std::fixed << val;
    }
    println(level, o.str());
}

void print_utils::vector(int level, const std::string &txt, const DoubleVector &val, int p, bool s) {
    if (p < 0) p = Printer::getPrecision();
    int w0 = Printer::getWidth() - 2;
    int w1 = w0 * 2 / 9;
    int w2 = w0 - 3 * w1;
    int w3 = w2 - (txt.size() + 1);

    std::stringstream o;
    o << " " << txt << std::string(w3, ' ') << ":";
    for (int i = 0; i < val.size(); i++) {
        if (std::isnan(val(i))) {
            o << std::setw(w1) << "N/A";
        } else if (s) {
            o << std::setw(w1) << std::setprecision(p) << std::scientific << val(i);
        } else {
            o << std::setw(w1) << std::setprecision(p) << std::fixed << val(i);
        }
    }
    println(level, o.str());
}

void print_utils::matrix(int level, const std::string &txt, const DoubleMatrix &val, int p, bool s) {
    if (p < 0) p = Printer::getPrecision();
    int w0 = Printer::getWidth() - 2;
    int w1 = w0 * 2 / 9;
    int w2 = w0 - 3 * w1;
    int w3 = w2 - (txt.size() + 1);

    std::stringstream o;
    for (int i = 0; i < val.rows(); i++) {
        if (i == 0) {
            o << " " << txt << std::string(w3, ' ') << ":";
        } else {
            o << " " << std::string(w2 - 1, ' ') << ":";
        }
        for (int j = 0; j < val.cols(); j++) {
            if (std::isnan(val(i, j))) {
                o << std::setw(w1) << "N/A";
            } else if (s) {
                o << std::setw(w1) << std::setprecision(p) << std::scientific << val(i, j);
            } else {
                o << std::setw(w1) << std::setprecision(p) << std::fixed << val(i, j);
            }
        }
        o << std::endl;
    }
    printout(level, o.str());
}

void print_utils::qmfunction(int level, const std::string &txt, const QMFunction &func, mrcpp::Timer &timer) {
    auto nodes = func.getNNodes(NUMBER::Total);
    auto memory = func.getSizeNodes(NUMBER::Total);
    auto time = timer.elapsed();
    mrcpp::print::tree(level, txt, nodes, memory, time);
}

} // namespace mrchem
