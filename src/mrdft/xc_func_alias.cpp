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
#include <xc_funcs.h>
#include <xc.h>

#include "xc_func_alias.h"

namespace mrdft {
    
void mapFunctionalName(std::string name, std::vector<int> &ids, std::vector<double> &coefs, double &customExx) {
    customExx = 0.0; // Allows for hard coded exact exchange
    // ensure name is upper case
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::toupper(c); });

    // LDA
    if (name == "SLATERX") {
        ids = {XC_LDA_X};
        coefs = {1.0};
        return;
    } else if (name == "VWN3C" || name == "VWN3") {
        ids = {XC_LDA_C_VWN_RPA};
        coefs = {1.0};
        return;
    } else if (name == "VWN5C" || name == "VWN5") {
        ids = {XC_LDA_C_VWN};
        coefs = {1.0};
        return;
    } else if (name == "SVWN3") {
        ids = {XC_LDA_C_VWN_RPA, XC_LDA_X};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "SVWN5" || name == "LDA") {
        ids = {XC_LDA_C_VWN, XC_LDA_X};
        coefs = {1.0, 1.0};
        return;

        // GGA
    } else if (name == "BECKEX") {
        ids = {XC_GGA_X_B88};
        coefs = {1.0};
        return;
    } else if (name == "BLYP") {
        ids = {XC_GGA_X_B88, XC_GGA_C_LYP};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "BP86") {
        ids = {XC_GGA_X_B88, XC_GGA_C_P86_FT};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "BPW91") {
        ids = {XC_GGA_X_B88, XC_GGA_C_PW91};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "OLYP") {
        ids = {XC_GGA_X_OPTX, XC_GGA_C_LYP};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "PBE") {
        ids = {XC_GGA_X_PBE, XC_GGA_C_PBE};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "KT1") {
        ids = {XC_GGA_XC_KT1};
        coefs = {1.0};
        return;
    } else if (name == "KT2") {
        //ok
        ids = {XC_GGA_XC_KT2};
        coefs = {1.0};
        return;
    } else if (name == "KT3") {
        //ok
        ids = {XC_GGA_XC_KT3};
        coefs = {1.0};
        return;
        
        // HYB GGA
    } else if (name == "B3LYP" || name == "B3LYP5") {
        ids = {XC_HYB_GGA_XC_B3LYP5};
        coefs = {1.0};
        return;
    } else if (name == "B3LYP-G") {
        //ok
        ids = {XC_HYB_GGA_XC_B3LYP};
        coefs = {1.0};
        return;
    } else if (name == "B3P86") {
        ids = {XC_HYB_GGA_XC_B3P86, XC_LDA_C_VWN_RPA, XC_LDA_C_VWN};
        coefs = {1.0, -1.0, 1.0};
        return;
    } else if (name == "B3P86-G") {
        //ok
        ids = {XC_HYB_GGA_XC_B3P86};
        coefs = {1.0};
        return;
    } else if (name == "PBE0") {
        ids = {XC_HYB_GGA_XC_PBEH};
        coefs = {1.0};
        return;

    // Other: Check if Libxc has this functional
    } else {
        int number = xc_functional_get_number(name.c_str());
        if (number == -1) { MSG_ABORT(name + " is not a known shorthand in MRChem nor a functional in Libxc.\n"); }

        ids = {number};
        coefs = {1.0};
        return;
    }
}

} // namespace mrdft
