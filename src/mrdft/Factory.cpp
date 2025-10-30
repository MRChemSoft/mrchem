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

#include "Factory.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <XCFun/xcfun.h>

#include "GGA.h"
#include "Grid.h"
#include "LDA.h"
#include "MRDFT.h"
#include "SpinGGA.h"
#include "SpinLDA.h"

namespace mrdft {

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA)
        : mra(MRA)
        , xcfun_p(xcfun_new(), xcfun_delete) {}

bool Factory::libxc;

// Todo: det finnes sikkert gode util funksjoner i libxc for dette
std::vector<int> Factory::mapFunctionalName(const std::string &name) const {
    // Map common functional names from XCFun to LibXC IDs, only LDAs for now
    if (name == "slaterx")                  return {XC_LDA_X};
    if (name == "svwn3")                    return {XC_LDA_X, XC_LDA_C_VWN_3};
    if (name == "svwn5c" || name == "VWN5") return {XC_LDA_C_VWN};
    if (name == "svwn5")                    return {XC_LDA_X, XC_LDA_C_VWN};
    if (name == "pbe")                      return {XC_GGA_X_PBE, XC_GGA_C_PBE};
    // if (name == "b3lyp-g")                  return {XC_HYB_GGA_XC_B3LYP3};
    if (name == "b3lyp-g")                  return {XC_LDA_X, XC_GGA_X_B88, XC_GGA_C_LYP, XC_LDA_C_VWN}; // not right
    if (name == "b3lyp")                    return {XC_HYB_GGA_XC_B3LYP};
    if (name == "pbe0")                     return {XC_HYB_GGA_XC_PBEH};
    if (name == "kt1")                      return {XC_LDA_X, XC_GGA_X_KT1, XC_LDA_C_VWN};

    std::cout << "!!!!! Add functional to mapFunctionalName(): " << name << std::endl;
    MSG_ABORT("Unknown functional for libxc")
    return {1};
}

void newMapFuncName(const std::string &name, std::vector<int> &ids, std::vector<double> &coeffs) {
    std::cout << "Name used in MapFunctionalName: " << name << std::endl;
    if (name == "pbe0") {
        // ids = {XC_GGA_X_PBE, XC_GGA_C_PBE}; // Both versions are the exact same
        // coeffs = {0.75, 1.0};
        ids = {XC_HYB_GGA_XC_PBEH};
        coeffs = {1.0};
        return;
    } else if (name == "slaterx" || name == "SLATERX") {
        ids = {XC_LDA_X};
        coeffs = {1.0};
        return;
    } else if (name == "BECKEX" || name == "beckex") {
        ids = {XC_GGA_X_B88};
        coeffs = {1.0};
        return;
    } else if (name == "svwn5c" || name == "VWN5C") {
        ids = {XC_LDA_C_VWN};
        coeffs = {1.0};
        return;
    } else if (name == "svwn5") {
        ids = {XC_LDA_C_VWN, XC_LDA_X};
        coeffs = {1.0, 1.0};
        return;
    } else if (name == "b3p86") {
        ids = {XC_HYB_GGA_XC_B3P86}; 
        coeffs = {1.0};
        return;
    } else if (name == "bpw91") {
        ids = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_C_PW91};
        coeffs = {.5, .5, 1.0}; // Closest
        return;
    } else {std::cout << "NO FUNC MAPPED" << std::endl;}
}

void Factory::setFunctional(const std::string &n, double c) {
    xcfun_set(xcfun_p.get(), n.c_str(), c);
    std::string name = n;
    std::cout << "xcfun func: " << n << std::endl;
    // std::vector<int> ids = this->mapFunctionalName(name);
    std::vector<int> ids;
    std::vector<double> coeffs;
    newMapFuncName(name, ids, coeffs);
    std::cout << "name sat" << std::endl;
    setLibxc(libxc);

    xc_func_type libxc_obj;
    for (size_t i = 0; i < ids.size(); i++) {
        if (spin) {
            if (xc_func_init(&libxc_obj, ids[i], XC_POLARIZED) != 0) {
                std::cout << "!!!!! Unknown functional (setfunctional)name : " << name << " id: " << ids[i] << "--" << xc_func_init(&libxc_obj, ids[i], XC_UNPOLARIZED) << std::endl;
            }
            xc_func_set_dens_threshold(&libxc_obj, cutoff);
            libxc_objects.push_back(libxc_obj);
            libxc_coeffs.push_back(c * coeffs[i]);
        } else {
            if (xc_func_init(&libxc_obj, ids[i], XC_UNPOLARIZED) != 0) {
                std::cout << "!!!!! Unknown functional (setfunctional)name : " << name << " id: " << ids[i] << "--" << xc_func_init(&libxc_obj, ids[i], XC_UNPOLARIZED) << std::endl;
            }
            xc_func_set_dens_threshold(&libxc_obj, cutoff);

            std::cout << "Functional number: " << libxc_objects.size() << ": " << n << std::endl;
            libxc_objects.push_back(libxc_obj);
            libxc_coeffs.push_back(c * coeffs[i]);
        }
    }
}

/** @brief Build a MRDFT object from the currently defined parameters */
std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);
    setLibxc(libxc);
    setFunctional("BECKEX", 1.0);

    // Init XCFun
    bool gga = xcfun_is_gga(xcfun_p.get());
    bool lda = not(gga);
    unsigned int mode = 1;                    //!< only partial derivative mode implemented
    unsigned int func_type = (gga) ? 1 : 0;   //!< only LDA and GGA supported for now
    unsigned int dens_type = 1 + spin;        //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int laplacian = 0;               //!< no laplacian
    unsigned int kinetic = 0;                 //!< no kinetic energy density
    unsigned int current = 0;                 //!< no current density
    unsigned int exp_derivative = not(gamma); //!< use gamma or explicit derivatives
    if (not(gga)) exp_derivative = 0;         //!< fall back to gamma-type derivatives if LDA
    xcfun_user_eval_setup(xcfun_p.get(), order, func_type, dens_type, mode, laplacian, kinetic, current, exp_derivative);
    // xcfun_eval_setup(xcfun_p.get(), XC_N_NX_NY_NZ, XC_PARTIAL_DERIVATIVES, 1);

    // Init MW derivative
    if (gga) {
        if (diff_s == "bspline") diff_p = std::make_unique<mrcpp::BSOperator<3>>(mra, 1);
        if (diff_s == "abgv_00") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
        if (diff_s == "abgv_55") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.5, 0.5);
    }

    // Init XC functional
    std::unique_ptr<Functional> func_p{nullptr};
    if (spin) {
        if (gga) func_p = std::make_unique<SpinGGA>(order, xcfun_p, diff_p);
        if (lda) func_p = std::make_unique<SpinLDA>(order, xcfun_p);
    } else {
        if (gga) func_p = std::make_unique<GGA>(order, xcfun_p, diff_p);
        if (lda) func_p = std::make_unique<LDA>(order, xcfun_p);
    }

    if (func_p == nullptr) MSG_ABORT("Invalid functional type");

    func_p->set_libxc_functional_object(libxc_objects, libxc_coeffs);

    diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(diff_p);
    func_p->setLogGradient(log_grad);
    // func_p->setDensityCutoff(cutoff);
    std::cout << "cutoff: " << cutoff << std::endl; // != cutoff set in input
    func_p->setDensityCutoff(1e-6);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);

    return mrdft_p;
}

} // namespace mrdft
