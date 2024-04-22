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

#include "GPESolver.h"
#include "chemistry/PhysicalConstants.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

PBESolver::PBESolver(const Permittivity &e,
                     const DHScreening &k,
                     const Density &rho_nuc,
                     std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                     int kain_hist,
                     int max_iter,
                     bool dyn_thrs,
                     SCRFDensityType density_type)
        : GPESolver(e, rho_nuc, P, D, kain_hist, max_iter, dyn_thrs, density_type)
        , kappa(k) {}

void PBESolver::computePBTerm(mrcpp::ComplexFunction &V_tot, const double salt_factor, mrcpp::ComplexFunction &pb_term) {
    // create a lambda function for the sinh(V) term and multiply it with kappa and salt factor to get the PB term
    auto sinh_f = [salt_factor](const double &V) { return (salt_factor / (4.0 * mrcpp::pi)) * std::sinh(V); };
    resetComplexFunction(pb_term);
    mrcpp::ComplexFunction sinhV;
    sinhV.alloc(NUMBER::Real);
    mrcpp::map(this->apply_prec / 100, sinhV.real(), V_tot.real(), sinh_f);

    mrcpp::cplxfunc::multiply(pb_term, sinhV, this->kappa, this->apply_prec);
}

void PBESolver::computeGamma(mrcpp::ComplexFunction &potential, mrcpp::ComplexFunction &out_gamma) {

    auto d_V = mrcpp::gradient(*derivative, potential.real()); // FunctionTreeVector
    resetComplexFunction(out_gamma);

    for (int d = 0; d < 3; d++) {
        auto C_pin = this->epsilon.getCavity_p();
        mrcpp::AnalyticFunction<3> d_cav(C_pin->getGradVector()[d]);
        mrcpp::ComplexFunction cplxfunc_prod;
        mrcpp::cplxfunc::multiply(cplxfunc_prod, get_func(d_V, d), d_cav, this->apply_prec, 1);
        // add result into out_gamma
        if (d == 0) {
            mrcpp::cplxfunc::deep_copy(out_gamma, cplxfunc_prod);
        } else {
            out_gamma.add(1.0, cplxfunc_prod);
        }
    }

    out_gamma.rescale(std::log((epsilon.getValueIn() / epsilon.getValueOut())) * (1.0 / (4.0 * mrcpp::pi)));
    mrcpp::clear(d_V, true);

    // add PB term
    mrcpp::ComplexFunction pb_term;
    auto salt_factor = 1.0; // placeholder for now, want to change it wrt to convergence in future
    computePBTerm(potential, salt_factor, pb_term);
    out_gamma.add(-1.0, pb_term);
}

} // namespace mrchem
