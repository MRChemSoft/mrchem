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

#include <memory>
#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/Printer>
#include <XCFun/xcfun.h>

#include "XCFun.h"

namespace mrdft {

XCFun::~XCFun() {
    if (xcfun != nullptr) {
        xcfun_delete(xcfun);
        xcfun = nullptr;
    }
}

void XCFun::setFunctional(const std::string &name, double c, double cutoff) {
    addFunctionalSpec(name, c);
    xcfun_set(xcfun, name.c_str(), c);
}

void XCFun::initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, int order, bool gamma) {
    auto *xcfun_p = xcfun;
    gga = xcfun_is_gga(xcfun_p);
    mgga = xcfun_is_metagga(xcfun_p);
    if (not(gga)) { if (not(mgga)) {lda = true; } }
    unsigned int mode = 1;                    //!< only partial derivative mode implemented
    unsigned int dens_type = 1 + spin;        //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int laplacian = 0;               //!< no laplacian
    unsigned int kinetic = mgga ? 1u : 0u;    //!< request kinetic energy density for meta-GGAs
    unsigned int current = 0;                 //!< no current density
    unsigned int exp_derivative = not(gamma); //!< use gamma or explicit derivatives
    if (lda) exp_derivative = 0;              //!< fall back to gamma-type derivatives if LDA
    unsigned int func_type = 0;               //!< LDA = 0, GGA = 1, mGGA = 2
        if (mgga) {
            func_type = 2;
        } else if (gga) {
            func_type = 1;
        } else {
            func_type = 0;
        }
    xcfun_user_eval_setup(xcfun_p, order, func_type, dens_type, mode, laplacian, kinetic, current, exp_derivative);
}

double XCFun::getAmountExx() const {
    double exx = 0.0;
    xcfun_get(xcfun, "exx", &exx);
    return exx;
}

void XCFun::printFunctionalReference(int out_txt_width) const {
    // Print header and provide wrapping utility via XClib helpers
    XCLib::printReferenceHeader();
    printout(0, xcfun_splash());
    std::cout << "\nXCFun functionals used in this calculation:\n";
    for (const auto &func_name : xcfun_func_names) {
        std::string xcfun_ref = xcfun_describe_long(func_name.c_str());
        std::string xcfun_ref_str = "  - " + xcfun_ref;
        XCLib::printWrap(xcfun_ref_str, out_txt_width, 4);
    }
    return;
}

void XCFun::callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts, int nInp, int nOut, double cutoff) const {
    if (nInp != xcfun_input_length(xcfun) or nOut != xcfun_output_length(xcfun)) { throw std::logic_error("Dimension mismatch!\n"); }

    for (int i = 0; i < nPts; i++) {
        if (spin) {
            if (inp(0, i) < cutoff and inp(1, i) < cutoff) continue;
        } else {
            if (inp(0, i) < cutoff) continue;
        }
        xcfun_eval(xcfun, inp.col(i).data(), out.col(i).data());
    }
}

} // mrdft
