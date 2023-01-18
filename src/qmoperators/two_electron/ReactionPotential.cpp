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
#include <MRCPP/Timer>

#include "ReactionPotential.h"

#include "qmfunctions/qmfunction_utils.h"

using SCRF_p = std::unique_ptr<mrchem::SCRF>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

ReactionPotential::ReactionPotential(SCRF_p scrf_p, OrbitalVector_p Phi_p)
        : QMPotential(1, false)
        , helper(std::move(scrf_p))
        , Phi(Phi_p) {}

void ReactionPotential::setup(double prec) {
    setApplyPrec(prec);
    auto thrs = this->helper->setConvergenceThreshold(prec);
    mrcpp::Timer timer;
    auto plevel = mrcpp::Printer::getPrintLevel();
    mrcpp::print::separator(3, '=');
    print_utils::centered_text(3, "Building Reaction operator");
    this->helper->printParameters();
    mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
    mrcpp::print::value(3, "Threshold", thrs, "(abs)", 5);
    mrcpp::print::separator(3, '-');
    auto potential = this->helper->setup(prec, this->Phi);
    qmfunction::deep_copy(*this, potential);
    if (plevel == 2) print_utils::qmfunction(2, "Reaction operator", *this, timer);
    mrcpp::print::footer(3, timer, 2);
}

void ReactionPotential::clear() {
    clearApplyPrec();
    this->helper->clear();
}

} // namespace mrchem
