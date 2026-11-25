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

#include "GGA.h"
#include "Grid.h"
#include "LDA.h"
#include "MRDFT.h"
#include "mGGA.h"

#ifndef DISABLE_LIBXC
#include "xc_libraries/Libxc.h"
#endif

#ifndef DISABLE_XCFUN
#include "xc_libraries/XCFun.h"
#endif

namespace mrdft {

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA, bool spin, const std::string &xclibname)
        : spin(spin), mra(MRA) {
    if (xclibname == "xcfun") {
#ifdef DISABLE_XCFUN
        MSG_ABORT("XCFun support disabled during compilation!");
#else
        xclib = std::make_unique<XCFun>(spin);
#endif
    } else if (xclibname == "libxc") {
#ifdef DISABLE_LIBXC
        MSG_ABORT("LibXC support disabled during compilation!");
#else
        xclib = std::make_unique<Libxc>(spin);
#endif
    } else if (xclibname == "auto") {
#ifndef DISABLE_XCFUN
        xclib = std::make_unique<XCFun>(spin);
#else
#ifndef DISABLE_LIBXC
        xclib = std::make_unique<Libxc>(spin);
#else
        MSG_ABORT("No XC library compiled in. DFT, and the sad/nao initial guesses (also for Hartree-Fock), need XCFun or Libxc; use guess_type gto or core (not possible with pseudopotentials, which always use the nao guess).");
#endif
#endif
    } else {
        MSG_ABORT("Unknown xc_library '" << xclibname << "' (use auto, xcfun or libxc)");
    }
}

std::string Factory::getXCLibName() const {
    return xclib->name();
}

void Factory::setFunctional(const std::string &name, double c) {
    xclib->setFunctional(name, c);
}

void Factory::setDensityCutoff(double c) {
    cutoff = c;
    xclib->setCutoff(c);
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
    // makepot() differentiates node by node (DerivativeCalculator::calcNode), which is only valid for a
    // bandwidth-0 operator. ABGV(0,0) is the only one.
    if (diff_s != "abgv_00") MSG_ABORT("Factory::build: the XC derivative must be abgv_00, got '" << diff_s << "'");
    diff_p = std::make_shared<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);

    // Init XC functional (meta-GGA first: a meta-GGA functional may also report gga)
    std::unique_ptr<Functional> func_p{nullptr};
    if (mgga) {
        if (spin) MSG_ABORT("Factory::build: spin-polarized meta-GGA is not supported yet"); // TODO: implement spin
        func_p = std::make_unique<mGGA>(order, spin, xclib);
    } else if (gga) {
        func_p = std::make_unique<GGA>(order, spin, xclib);
    } else if (lda) {
        func_p = std::make_unique<LDA>(order, spin, xclib);
    } else {
        MSG_ABORT("Case not handled");
    }
    if (func_p == nullptr) MSG_ABORT("Invalid functional type");
    func_p->setDerivOp(diff_p);
    if (mgga) {
        // TODO: make optional derivative operator from input preference (eg xc_derivOp/tau/...)
        // tau, the orbital grad rho and the tau operator need an antisymmetric derivative (bspline, abgv_55, ph(?))
        func_p->setTauDerivOp(std::make_shared<mrcpp::BSOperator<3>>(mra, 1));
    }
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft
