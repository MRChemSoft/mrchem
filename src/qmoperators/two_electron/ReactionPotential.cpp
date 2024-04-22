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
#include "qmfunctions/density_utils.h"

using GPESolver_p = std::unique_ptr<mrchem::GPESolver>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

ReactionPotential::ReactionPotential(GPESolver_p gpesolver_p, OrbitalVector_p Phi, bool mpi_share)
        : QMPotential(1, mpi_share)
        , solver(std::move(gpesolver_p))
        , orbitals(Phi) {}

void ReactionPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    auto thrs = this->solver->setConvergenceThreshold(prec);
    mrcpp::Timer timer;
    auto plevel = mrcpp::Printer::getPrintLevel();
    mrcpp::print::separator(3, '=');
    print_utils::centered_text(3, "Building Reaction operator");
    this->solver->printParameters();
    mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
    mrcpp::print::value(3, "Threshold", thrs, "(abs)", 5);
    mrcpp::print::separator(3, '-');
    auto potential = this->computePotential(prec);
    mrcpp::cplxfunc::deep_copy(*this, potential);
    if (plevel == 2) print_utils::qmfunction(2, "Reaction operator", *this, timer);
    mrcpp::print::footer(3, timer, 2);
}

void ReactionPotential::clear() {
    mrcpp::ComplexFunction::free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();
    this->solver->clear();
}

} // namespace mrchem
