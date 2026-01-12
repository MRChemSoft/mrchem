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

void MapFuncName(std::string name, std::vector<int> &ids, std::vector<double> &coefs) {
    // ensure name is upper case
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::toupper(c); });

    std::cout << "Name used in MapFunctionalName: " << name << std::endl;
    if (name == "PBE0") {
        ids = {XC_HYB_GGA_XC_PBEH};
        coefs = {1.0};
        return;
    } else if (name == "SLATERX") {
        ids = {XC_LDA_X};
        coefs = {1.0};
        return;
    } else if (name == "BECKEX") {
        ids = {XC_GGA_X_B88};
        coefs = {1.0};
        return;
    } else if (name == "VWN5C") {
        ids = {XC_LDA_C_VWN};
        coefs = {1.0};
        return;
    } else if (name == "SVWN5") {
        ids = {XC_LDA_C_VWN, XC_LDA_X};
        coefs = {1.0, 1.0};
        return;
    } else if (name == "B3P86") {
        ids = {XC_HYB_GGA_XC_B3P86};
        coefs = {1.0};
        return;
    } else if (name == "BPW91") {
        ids = {XC_GGA_X_B88, XC_GGA_C_PW91};
        coefs = {1.0, 1.0};
        return;
    } else {
        // Change any dashes to underscores
        for (size_t i = 0; i < name.size(); i++) {
            if (name[i] == '-') { name[i] = '_'; }
        }

        // Check if Libxc has this functional
        int number = xc_functional_get_number(name.c_str());
        if (number == -1) { throw std::logic_error("Got name " + name + " but this is not a known shorthand in MRChem nor a functional in Libxc\n"); }

        ids = {number};
        coefs = {1.0};
        return;
    }
}

void Factory::setFunctional(const std::string &name, double c) {
    setLibxc(libxc); // should probably be where setFunctional is called

    if (not libxc) {
        xcfun_set(xcfun_p.get(), name.c_str(), c);

    } else {
        std::vector<int> ids;
        std::vector<double> coefs;

        MapFuncName(name, ids, coefs);
        xc_func_type libxc_obj;
        for (size_t i = 0; i < ids.size(); i++) {
            auto return_code = xc_func_init(&libxc_obj, ids[i], spin ? XC_POLARIZED : XC_UNPOLARIZED);
            if (return_code != 0) { std::cout << "!!!!! Unknown functional (setfunctional)name : " << name << " id: " << ids[i] << "--" << return_code << std::endl; }
            xc_func_set_dens_threshold(&libxc_obj, cutoff);

            std::cout << "Functional number: " << libxc_objects.size() << ": " << name << std::endl;
            libxc_objects.push_back(libxc_obj);
            libxc_coefs.push_back(c * coefs[i]);
        }
    }
}

/** @brief Build a MRDFT object from the currently defined parameters */
std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);
    setLibxc(libxc);

    // Init XCFun or Libxc
    bool gga;
    if (not libxc) {
        gga = xcfun_is_gga(xcfun_p.get());
        unsigned int mode = 1;                    //!< only partial derivative mode implemented
        unsigned int func_type = (gga) ? 1 : 0;   //!< only LDA and GGA supported for now
        unsigned int dens_type = 1 + spin;        //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
        unsigned int laplacian = 0;               //!< no laplacian
        unsigned int kinetic = 0;                 //!< no kinetic energy density
        unsigned int current = 0;                 //!< no current density
        unsigned int exp_derivative = not(gamma); //!< use gamma or explicit derivatives
        if (not(gga)) exp_derivative = 0;         //!< fall back to gamma-type derivatives if LDA
        xcfun_user_eval_setup(xcfun_p.get(), order, func_type, dens_type, mode, laplacian, kinetic, current, exp_derivative);

    } else {
        for (const auto &f : libxc_objects) switch (f.info->family) {
                case XC_FAMILY_LDA:
#ifdef XC_FAMILY_HYB_GGA
                case XC_FAMILY_HYB_LDA:
#endif
                    gga = false;
                    break;

                case XC_FAMILY_GGA:
#ifdef XC_FAMILY_HYB_GGA
                case XC_FAMILY_HYB_GGA:
#endif
                    gga = true;
                    break;

                case XC_FAMILY_MGGA:
                case XC_FAMILY_HYB_MGGA:
                    gga = false; // eliminate unused variable warning
                    MSG_ABORT("Meta-GGA functionals are not supported in MRChem.!\n");

                default:
                    gga = false; // eliminate unused variable warning
                    MSG_ABORT("Case not handled in MRChem!\n");
            }
    }
    bool lda = not gga;

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
    if (libxc) { func_p->set_libxc_functional_object(libxc_objects, libxc_coefs); }
    diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(diff_p);
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft
