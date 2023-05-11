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

/** @brief Build a MRDFT object from the currently defined parameters */
std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);

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
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft
