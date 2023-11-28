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

#include "LPBESolver.h"

#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "GPESolver.h"
#include "chemistry/PhysicalConstants.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmoperators/two_electron/ReactionPotential.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {
LPBESolver::LPBESolver(const Permittivity &e,
                       const DHScreening &k,
                       const Density &rho_nuc,
                       std::shared_ptr<mrcpp::PoissonOperator> P,
                       std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                       int kain_hist,
                       int max_iter,
                       bool dyn_thrs,
                       SCRFDensityType density_type)
        : PBESolver(e, k, rho_nuc, P, D, kain_hist, max_iter, dyn_thrs, density_type) {}
// TODO separate this for the linear and non-linear solver
void LPBESolver::computePBTerm(mrcpp::ComplexFunction &V_tot, const double salt_factor, mrcpp::ComplexFunction &pb_term) {
    resetComplexFunction(pb_term);
    mrcpp::cplxfunc::multiply(pb_term, V_tot, this->kappa, this->apply_prec);
    pb_term.rescale(salt_factor / (4.0 * mrcpp::pi));
}

} // namespace mrchem
