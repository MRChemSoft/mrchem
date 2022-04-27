/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <locale>

#include <MRCPP/Printer>
#include <MRCPP/Timer>

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

void print_utils::json(int level, const nlohmann::json &j, bool ralign) {
    // Determine longest name
    // This will be used for the spacing left of :
    // if the name is too long to be aligned with
    // the other sections
    int lshift = 0;
    for (const auto &item : j.items()) {
        if (item.key().size() > lshift) lshift = item.key().size();
    }

    // Loop over json items
    for (const auto &item : j.items()) {
        // Extract key and value from json
        std::string key = item.key();
        std::stringstream o_val;
        o_val << item.value();
        std::string val = o_val.str();

        // Remove quotes from val
        val.erase(std::remove(val.begin(), val.end(), '\"'), val.end());

        // Standard shift to align all colons
        int w0 = Printer::getWidth() - 2; // Defines the printable width
        int w1 = w0 * 2 / 9;              // Space dedicated to the json key
        int w2 = w0 - 3 * w1;
        int w3 = w2 - (val.size() + 1);

        // Some paddings for book keeping
        int frontEndPadding = 2; // Empty space at beginning and end
        int colonPadding = 3;    // Two empty spaces around a single colon

        // Use standard spacing if longest name fits
        if (w2 > lshift) lshift = w2 - frontEndPadding;

        // Calculate the shift needed for right-aligning
        int rshift = (ralign) ? Printer::getWidth() - lshift - val.size() - frontEndPadding - colonPadding : 0;

        // Check that rshift is not negative (caused by very long names)
        if (rshift < 0) rshift = 0;

        // Print line
        std::printf(" %-*s%-s%-s%-s \n", lshift, key.c_str(), " : ", std::string(rshift, ' ').c_str(), val.c_str());
    }
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

// trim from start (in place)
void print_utils::ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
void print_utils::rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

// trim from both ends (in place)
void print_utils::trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
std::string print_utils::ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
std::string print_utils::rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
std::string print_utils::trim_copy(std::string s) {
    trim(s);
    return s;
}

} // namespace mrchem
