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

void mapFunctionalName(std::string name, std::vector<int> &ids, std::vector<double> &coefs) {
    // ensure name is upper case
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::toupper(c); });

    // functionals defined in mrchem not explicitly translated to libxc:
    // b3lyp-g     {"b3lyp-g", "Becke-3-paramater-LYP (VWN3 form)", {{"slaterx", 0.80}, {"beckecorrx", 0.72}, {"lypc", 0.81}, {"vwn3c", 0.19}, {"exx", 0.20}}},
    // b3p86-g     {"b3p86-g", "Becke-3-paramater-LYP (VWN3 form)", {{"slaterx", 0.80}, {"beckecorrx", 0.72}, {"p86corrc", 0.81}, {"vwn3c", 1.0}, {"exx", 0.20}}}

    // LDA
    if (name == "SLATERX") {
        ids = {XC_LDA_X};
        coefs = {1.0};
        return;
    } else if (name == "VWN3C" || name == "VWN3") {
        ids = {XC_LDA_C_VWN_3};
        coefs = {1.0};
        return;
    } else if (name == "VWN5C" || name == "VWN5") {
        ids = {XC_LDA_C_VWN};
        coefs = {1.0};
        return;
    } else if (name == "SVWN3") {
        ids = {XC_LDA_C_VWN_3, XC_LDA_X};
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
        // xcfun def:     {"blyp", "Becke exchange and LYP correlation", {{"beckex", 1.0}, {"lypc", 1.0}}}
        ids = {XC_GGA_X_B88, XC_GGA_C_LYP};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "BP86") {
        // xcfun def:     {"bp86", "Becke-Perdew 1986", {{"beckex", 1.0}, {"p86c", 1.0}}}
        ids = {XC_GGA_X_B88, XC_GGA_C_P86};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "BPW91") {
        // xcfun def:     {"bpw91", "Becke 88 exchange+PW91", {{"beckex", 1.0}, {"pw91c", 1.0}}}
        ids = {XC_GGA_X_B88, XC_GGA_C_PW91};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "OLYP") {
        // xcfun def:     {"olyp", "LYP correlation and OPTX exchange", {{"lypc", 1.0}, {"optx", 1.0}}}
        ids = {XC_GGA_X_OPTX, XC_GGA_C_LYP};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "PBE") {
        // NB: not the exact same parameters, eq to 1e-7 for H2
        ids = {XC_GGA_X_PBE, XC_GGA_C_PBE};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "KT1") {
        // xcfun def:     {"kt1", "Keal-Tozer 2", {{"slaterx", 1.}, {"ktx", -0.006}, {"vwn5c", 1.}}}
        // libxc def:   static int   funcs_id  [2] = {XC_GGA_X_KT1, XC_LDA_C_VWN}; static double funcs_coef[2] = {1.0, 1.0};
        // if 1 * slaterx - 0.006 ktx = 1 XC_GGA_X_KT1 -> EQUAL
        ids = {XC_GGA_XC_KT1};
        coefs = {1.0};
        return;
    } else if (name == "KT2") {
        // xcfun def:     {"kt2", "Keal-Tozer 2", {{"slaterx", 1.07173}, {"ktx", -0.006}, {"vwn5c", 0.576727}}}
        // libxc def:   static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_KT1, XC_LDA_C_VWN}; static double funcs_coef[3] = {1.07173 - 1.0, 1.0, 0.576727}
        ids = {XC_GGA_XC_KT2};
        coefs = {1.0};
        return;
    } else if (name == "KT3") {
        // xcfun def :     {"kt3", "Keal-Tozer 3", {{"slaterx", 1.092}, {"ktx", -0.004}, {"optxcorr", -0.925452}, {"lypc", 0.864409}}}
        // libxc def:   static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_C_LYP, XC_GGA_X_KT1, XC_GGA_X_OPTX}; double funcs_coef[4] = { alpha - eps*a_optx/b_optx - 1.0, beta, 1, eps/b_optx};  
        ids = {XC_GGA_XC_KT3};
        coefs = {1.0};
        return;

    // HYB GGA
    } else if (name == "B3LYP") {
        // Keep as b3lyp5 for now to be consistent with xcfun
        // TODO: change the definition of b3lyp in mrchem to not be b3lyp5
        ids = {XC_HYB_GGA_XC_B3LYP5};
        // ids = {XC_HYB_GGA_XC_B3LYP};
        coefs = {1.0};
        return;
    } else if (name == "B3LYP5") {
        // libxc def:   static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN, XC_GGA_C_LYP}; funcs_coefs = set by ext_param
        ids = {XC_HYB_GGA_XC_B3LYP5};
        coefs = {1.0};
        return;
    } else if (name == "B3P86") {
        // NB: not the exact same parameters, eq to 1e-3 for H2
        // xcfun def:     {"b3p86", "Becke-3-paramater-LYP (VWN5 form)", {{"slaterx", 0.80}, {"beckecorrx", 0.72}, {"p86corrc", 0.81}, {"vwn5c", 1.0}, {"exx", 0.20}}}
        // libxc def:   static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_P86}; funcs_coefs = set by ext_param
        // different vwn: vwn5 vs vwn rpa
        // change to be similar to xcfun?
        ids = {XC_HYB_GGA_XC_B3P86};
        coefs = {1.0};
        return;
    } else if (name == "PBE0") {
        // NB: not the exact same parameters, equivalent to 1e-7 for H2
        ids = {XC_HYB_GGA_XC_PBEH};
        coefs = {1.0};
        return;

    // Other: Check if Libxc has this functional
    } else {
        int number = xc_functional_get_number(name.c_str());
        if (number == -1) { MSG_ABORT(name + " is not a known shorthand in MRChem nor a functional in Libxc!\n"); }

        ids = {number};
        coefs = {1.0};
        return;
    }
}

} // namespace mrdft
