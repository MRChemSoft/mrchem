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

#include "PBESolver.h"

#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "chemistry/PhysicalConstants.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmoperators/two_electron/ReactionPotential.h"
#include "scf_solver/KAIN.h"
#include "utils/print_utils.h"

#include <nlohmann/json.hpp>

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

PBESolver::PBESolver(Permittivity e,
                     DHScreening k,
                     const Nuclei &N,
                     PoissonOperator_p P,
                     DerivativeOperator_p D,
                     double orb_prec,
                     int kain_hist,
                     int max_iter,
                     bool acc_pot,
                     bool dyn_thrs,
                     std::string density_type)
        : GPESolver(e, N, P, D, orb_prec, kain_hist, max_iter, acc_pot, dyn_thrs, density_type)
        , do_static_salt(false)
        , kappa(false) {
    mrcpp::cplxfunc::project(this->kappa, k, NUMBER::Real, this->apply_prec);
}

void PBESolver::computeGamma(mrcpp::ComplexFunction &potential, mrcpp::ComplexFunction &out_gamma) {
    auto d_V = mrcpp::gradient(*derivative, potential.real());
    resetComplexFunction(out_gamma);
    mrcpp::dot(this->apply_prec, out_gamma.real(), d_V, this->d_cavity);
    out_gamma.rescale(std::log((epsilon.getEpsIn() / epsilon.getEpsOut())) * (1.0 / (4.0 * mrcpp::pi)));
    mrcpp::clear(d_V, true);
    computePBTerm(potential, 1.0);
    mrcpp::cplxfunc::add(out_gamma, 1.0, out_gamma, -1.0, this->pbe_term, -1.0); // this adds the PB term to the gamma
}

// TODO separate this for the linear and non-linear solver
void PBESolver::computePBTerm(mrcpp::ComplexFunction &V_tot,
                              double salt_factor) { // make a lambda function which evaluates std::sinh(V_tot) and multiplies it with this->kappa for the poisson-boltzmann equation
    auto sinh_f = [salt_factor](const double &V) { return (salt_factor / (4.0 * mrcpp::pi)) * std::sinh(V); };
    resetComplexFunction(this->pbe_term);
    mrcpp::ComplexFunction sinhV;
    sinhV.alloc(NUMBER::Real);
    mrcpp::map(this->apply_prec / 100, sinhV.real(), V_tot.real(), sinh_f);

    mrcpp::cplxfunc::multiply(this->pbe_term, this->kappa, sinhV, this->apply_prec);
}

} // namespace mrchem
