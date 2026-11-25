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

#include "Libxc.h"
#include <MRCPP/Printer>
namespace mrdft {

void Libxc::mapFunctionalName(std::string name, std::vector<int> &ids, std::vector<double> &coefs, double &customExx) {
    customExx = 0.0; // Allows for hard coded exact exchange
    // ensure name is upper case
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::toupper(c); });

    // LDA
    if (name == "SLATERX" || name == "SLATER" || name == "LDAX") {
        ids = {XC_LDA_X};
        coefs = {1.0};
        return;
    } else if (name == "VWN3C" || name == "VWN3") {
        ids = {XC_LDA_C_VWN_RPA};
        coefs = {1.0};
        return;
    } else if (name == "VWN5C" || name == "VWN5" || name == "VWN") {
        ids = {XC_LDA_C_VWN};
        coefs = {1.0};
        return;
    } else if (name == "SVWN3") {
        ids = {XC_LDA_C_VWN_RPA, XC_LDA_X};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "SVWN5" || name == "SVWN" || name == "LDA") {
        ids = {XC_LDA_C_VWN, XC_LDA_X};
        coefs = {1.0, 1.0};
        return;

        // GGA
    } else if (name == "BECKEX" || name == "B88X") {
        ids = {XC_GGA_X_B88};
        coefs = {1.0};
        return;
    } else if (name == "PBEX") {
        ids = {XC_GGA_X_PBE};
        coefs = {1.0};
        return;
    } else if (name == "OPTX") {
        ids = {XC_GGA_X_OPTX};
        coefs = {1.0};
        return;
    } else if (name == "LYP") {
        ids = {XC_GGA_C_LYP};
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
    } else if (name == "B3P86-XCFUN") {
        // Same as XCFun: 0.08 Slater + 0.72 B88 + 1.0 VWN5 + 0.81 (P86 - PZ81) + 0.20 exx
        // (XC_HYB_GGA_XC_B3P86 has 0.19*VWN_RPA + 0.81*P86, i.e. PZ81 local correlation;
        //  P86_FT is used, as in BP86, since it matches XCFun p86c)
        ids = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN, XC_GGA_C_P86_FT, XC_LDA_C_PZ};
        coefs = {0.08, 0.72, 1.0, 0.81, -0.81};
        customExx = 0.20;
        return;
    } else if (name == "B3P86-G-XCFUN") {
        // Same as XCFun: B3P86-XCFUN but with VWN3 (VWN_RPA) local correlation
        ids = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_P86_FT, XC_LDA_C_PZ};
        coefs = {0.08, 0.72, 1.0, 0.81, -0.81};
        customExx = 0.20;
        return;
    } else if (name == "PBE0" || name == "PBEH") {
        ids = {XC_HYB_GGA_XC_PBEH};
        coefs = {1.0};
        return;
    } else if (name == "B97") {
        ids = {XC_HYB_GGA_XC_B97};
        coefs = {1.0};
        return;
    } else if (name == "B97-1" || name == "B97_1") {
        ids = {XC_HYB_GGA_XC_B97_1};
        coefs = {1.0};
        return;
    } else if (name == "B97-2" || name == "B97_2") {
        ids = {XC_HYB_GGA_XC_B97_2};
        coefs = {1.0};
        return;

        // META-GGA
    } else if (name == "TPSS") {
        ids = {XC_MGGA_X_TPSS, XC_MGGA_C_TPSS};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "REVTPSS") {
        ids = {XC_MGGA_X_REVTPSS, XC_MGGA_C_REVTPSS};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "M06-L" || name == "M06L") {
        ids = {XC_MGGA_X_M06_L, XC_MGGA_C_M06_L};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "SCAN") {
        ids = {XC_MGGA_X_SCAN, XC_MGGA_C_SCAN};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "RSCAN") {
        ids = {XC_MGGA_X_RSCAN, XC_MGGA_C_RSCAN};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "RPPSCAN" || name == "R++SCAN") {
        ids = {XC_MGGA_X_RPPSCAN, XC_MGGA_C_RPPSCAN};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "R2SCAN") {
        ids = {XC_MGGA_X_R2SCAN, XC_MGGA_C_R2SCAN};
        coefs = {1.0, 1.0};
        return;

        // HYB META-GGA
    } else if (name == "TPSSH") {
        ids = {XC_HYB_MGGA_XC_TPSSH};
        coefs = {1.0};
        return;
    } else if (name == "REVTPSSH") {
        ids = {XC_HYB_MGGA_XC_REVTPSSH};
        coefs = {1.0};
        return;
    } else if (name == "TPSS0") {
        ids = {XC_HYB_MGGA_XC_TPSS0};
        coefs = {1.0};
        return;
    } else if (name == "M06") {
        ids = {XC_HYB_MGGA_X_M06, XC_MGGA_C_M06};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "M06-2X" || name == "M062X") {
        ids = {XC_HYB_MGGA_X_M06_2X, XC_MGGA_C_M06_2X};
        coefs = {1.0, 1.0};
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
