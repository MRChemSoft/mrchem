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

#include "MRDFT.h"
#include "Grid.h"
#include "LDA.h"
#include "GGA.h"
#include "mGGA.h"
#include "SpinGGA.h"
#include "SpinLDA.h"
#include "xc_libraries/Libxc.h"
#include "xc_libraries/XCFun.h"

namespace mrdft {

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA, bool spin_enabled, const std::string &xclibname)
    : spin_enabled(spin_enabled), xclibname(xclibname), mra(MRA) {
    if (xclibname == "xcfun") {
        xclib = std::make_unique<XCFun>(spin_enabled);
    } else if (xclibname == "libxc") {
        xclib = std::make_unique<Libxc>(spin_enabled);
    }
}

void Factory::setFunctional(const std::string &name, double c) {
    xclib->setFunctional(name, c, cutoff);
}

std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);

    // Init Libxc or XCFun and set functional family bools
    bool gga = false;
    bool lda = false;
    bool mgga = false;
    xclib->initFunctionalLibrary(lda, gga, mgga, order, gamma);

    // Init MW derivative
    if (gga || mgga) {
        if (!diff_p) {
            if (diff_s == "bspline") {
                diff_p = std::make_shared<mrcpp::BSOperator<3>>(mra, 1);
            } else if (diff_s == "abgv_00") {
                diff_p = std::make_shared<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
            } else if (diff_s == "abgv_55") {
                diff_p = std::make_shared<mrcpp::ABGVOperator<3>>(mra, 0.5, 0.5);
            } else {
                MSG_ABORT("Factory::build: unknown derivative operator '" + diff_s + "'");
            }
        }
    }

    std::unique_ptr<mrchem::KineticOperator> kin_p{nullptr};
    if (mgga) {
        if (spin_enabled) MSG_ABORT("Factory::build: spin-polarized meta-GGA is not supported yet");
        kin_p = std::make_unique<mrchem::KineticOperator>(diff_p);
    }

    std::unique_ptr<Functional> func_p{nullptr};
    if (spin_enabled) {
        if (mgga) {
            MSG_ABORT("Factory::build: spin-polarized meta-GGA is not supported yet");
        } else if (gga) {
            func_p = std::make_unique<SpinGGA>(order, xclib, diff_p);
        } else if (lda) {
            func_p = std::make_unique<SpinLDA>(order, xclib);
        } else {
            MSG_ABORT("Case not handled");
        }
    } else {
        if (mgga) {
            func_p = std::make_unique<mGGA>(order, xclib, diff_p);
        } else if (gga) {
            func_p = std::make_unique<GGA>(order, xclib, diff_p);
        } else if (lda) {
            func_p = std::make_unique<LDA>(order, xclib);
        } else {
            MSG_ABORT("Case not handled");
        }
    }

    if (func_p == nullptr) MSG_ABORT("Invalid functional type");
    diff_p = std::make_shared<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(diff_p);
    func_p->setKinOp(kin_p);
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft
