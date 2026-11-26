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

#include <MRCPP/Printer>
#include "XCLib.h"

namespace mrdft {


void XCLib::printReferenceHeader() {
    // Only run on main thread
    if (mrcpp::mpi::world_rank != 0) {
        return;
    }

    auto pwidth = mrcpp::Printer::getWidth();
    int txt_width = 50;
    auto pre_spaces = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');

    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '*');
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << "                  XC Functional                   " << post_str);
    println(0, pre_str << "                                                  " << post_str);
    mrcpp::print::separator(0, '*', 1);
}

void XCLib::printWrap(const std::string &str_in, std::size_t txt_width, int indent) {
    std::string str = str_in;
    const std::string continuation_indent(indent, ' ');
    size_t offset = 0;
    while (offset + txt_width < str.size()) {
        size_t newline_pos = str.find('\n', offset);
        if (newline_pos < offset + txt_width) {
            if (newline_pos != std::string::npos && newline_pos + 1 < str.size()) {
                str.insert(newline_pos + 1, continuation_indent);
                offset = newline_pos + 1 + continuation_indent.size();
            } else {
                offset = newline_pos + 1;
            }
            continue;
        }

        size_t space_pos = str.rfind(' ', offset + txt_width);
        if (space_pos == std::string::npos || space_pos - offset > txt_width) {
            space_pos = offset + txt_width;
            str.insert(space_pos, "\n" + continuation_indent);
        } else {
            str[space_pos] = '\n';
            str.insert(space_pos + 1, continuation_indent);
        }
        offset = space_pos + 1 + continuation_indent.size();
    }
    std::cout << str;
}

} // mrdft
