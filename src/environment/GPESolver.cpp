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

#include "GPESolver.h"

#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "Permittivity.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmoperators/two_electron/ReactionPotential.h"
#include "scf_solver/KAIN.h"
#include "utils/print_utils.h"

#include <nlohmann/json.hpp>

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;

namespace mrchem {

GPESolver::GPESolver(const Permittivity &e, const Density &rho_nuc, PoissonOperator_p P, DerivativeOperator_p D, int kain_hist, int max_iter, bool dyn_thrs, SCRFDensityType density_type)
        : dynamic_thrs(dyn_thrs)
        , density_type(density_type)
        , max_iter(max_iter)
        , history(kain_hist)
        , epsilon(e)
        , rho_nuc(rho_nuc)
        , Vr_n(false)
        , derivative(D)
        , poisson(P) {}

GPESolver::~GPESolver() {
    this->rho_nuc.free(NUMBER::Real);
    clear();
}

void GPESolver::clear() {
    this->apply_prec = -1.0;
}

double GPESolver::setConvergenceThreshold(double prec) {
    // converge_thrs should be in the interval [prec, 1.0]
    this->conv_thrs = prec;
    if (this->dynamic_thrs and this->mo_residual > 10 * prec) this->conv_thrs = std::min(1.0, this->mo_residual);
    return this->conv_thrs;
}

void GPESolver::computeDensities(const Density &rho_el, Density &rho_out) {
    Timer timer;

    switch (this->density_type) {
        case SCRFDensityType::TOTAL:
            mrcpp::cplxfunc::add(rho_out, 1.0, rho_el, 1.0, this->rho_nuc, -1.0);
            break;
        case SCRFDensityType::ELECTRONIC:
            // using const_cast is absolutely EVIL here and in principle shouldn't be needed!
            // however MRCPP is not everywhere const-correct, so here we go!
            mrcpp::cplxfunc::deep_copy(rho_out, const_cast<Density &>(rho_el));
            break;
        case SCRFDensityType::NUCLEAR:
            mrcpp::cplxfunc::deep_copy(rho_out, this->rho_nuc);
            break;
        default:
            MSG_ABORT("Invalid density type");
    }
    print_utils::qmfunction(3, "Vacuum density", rho_out, timer);
}

void GPESolver::computeGamma(mrcpp::ComplexFunction &potential, mrcpp::ComplexFunction &out_gamma) {

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
}

mrcpp::ComplexFunction GPESolver::solvePoissonEquation(const mrcpp::ComplexFunction &in_gamma, const Density &rho_el) {
    mrcpp::ComplexFunction Poisson_func;
    mrcpp::ComplexFunction rho_eff;
    mrcpp::ComplexFunction first_term;
    mrcpp::ComplexFunction Vr_np1;
    Vr_np1.alloc(NUMBER::Real);

    auto eps_inv_func = mrcpp::AnalyticFunction<3>([this](const mrcpp::Coord<3> &r) { return 1.0 / this->epsilon.evalf(r); });
    Density rho_tot(false);
    computeDensities(rho_el, rho_tot);

    mrcpp::cplxfunc::multiply(first_term, rho_tot, eps_inv_func, this->apply_prec);

    mrcpp::cplxfunc::add(rho_eff, 1.0, first_term, -1.0, rho_tot, -1.0);
    rho_tot.free(NUMBER::Real);

    mrcpp::cplxfunc::add(Poisson_func, 1.0, in_gamma, 1.0, rho_eff, -1.0);
    mrcpp::apply(this->apply_prec, Vr_np1.real(), *poisson, Poisson_func.real());
    return Vr_np1;
}

void GPESolver::accelerateConvergence(mrcpp::ComplexFunction &dfunc, mrcpp::ComplexFunction &func, KAIN &kain) {
    OrbitalVector phi_n(0);
    OrbitalVector dPhi_n(0);
    phi_n.push_back(Orbital(SPIN::Paired));
    dPhi_n.push_back(Orbital(SPIN::Paired));

    phi_n[0] = func;
    dPhi_n[0] = dfunc;

    phi_n[0].setRank(mrcpp::mpi::wrk_rank);
    dPhi_n[0].setRank(mrcpp::mpi::wrk_rank);

    kain.accelerate(this->apply_prec, phi_n, dPhi_n);

    func = phi_n[0];
    dfunc = dPhi_n[0]; // see if you can remove these two lines

    phi_n.clear();
    dPhi_n.clear();
}

void GPESolver::runMicroIterations(const mrcpp::ComplexFunction &V_vac, const Density &rho_el) {
    KAIN kain(this->history);
    kain.setLocalPrintLevel(10);

    mrcpp::print::separator(3, '-');
    auto update = 10.0, norm = 1.0;

    auto iter = 1;
    while (update >= this->conv_thrs && iter <= max_iter) {
        Timer t_iter;
        // solve the poisson equation
        mrcpp::ComplexFunction V_tot;
        mrcpp::ComplexFunction gamma_n;
        mrcpp::ComplexFunction dVr_n;

        mrcpp::cplxfunc::add(V_tot, 1.0, this->Vr_n, 1.0, V_vac, -1.0);
        computeGamma(V_tot, gamma_n);
        auto Vr_np1 = solvePoissonEquation(gamma_n, rho_el);
        norm = Vr_np1.norm();

        // use a convergence accelerator

        mrcpp::cplxfunc::add(dVr_n, 1.0, Vr_np1, -1.0, this->Vr_n, -1.0);
        update = dVr_n.norm();

        if (iter > 1 and this->history > 0) {
            accelerateConvergence(dVr_n, Vr_n, kain);
            Vr_np1.free(NUMBER::Real);
            mrcpp::cplxfunc::add(Vr_np1, 1.0, Vr_n, 1.0, dVr_n, -1.0);
        }

        // set up for next iteration
        resetComplexFunction(this->Vr_n);
        mrcpp::cplxfunc::deep_copy(this->Vr_n, Vr_np1);
        Vr_np1.free(NUMBER::Real);

        printConvergenceRow(iter, norm, update, t_iter.elapsed());

        // compute new PB term
        if (update < this->conv_thrs) break;
        iter++;
    }

    if (iter > max_iter) println(0, "Reaction potential failed to converge after " << iter - 1 << " iterations, residual " << update);
    mrcpp::print::separator(3, '-');

    kain.clear();
}

void GPESolver::printConvergenceRow(int i, double norm, double update, double time) const {
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 9;
    auto w2 = 2 * w0 / 9;
    auto w3 = w0 - w1 - 2 * w2;

    std::string time_unit = " sec";
    if (time < 0.01) {
        time = time * 1000.0;
        time_unit = "  ms";
    } else if (time > 60.0) {
        time = time / 60.0;
        time_unit = " min";
    }

    std::stringstream o_txt;
    o_txt << " Iter " << std::setw(w1 - 6) << i;
    o_txt << " : " << std::setw(w3 - 3) << std::setprecision(2 * pprec) << std::scientific << norm;
    o_txt << std::setw(w2) << std::setprecision(pprec) << std::scientific << update;
    o_txt << std::setw(w2 - 4) << std::setprecision(2) << std::fixed << time << time_unit;
    println(3, o_txt.str());
}

mrcpp::ComplexFunction &GPESolver::solveEquation(double prec, const Density &rho_el) {
    this->apply_prec = prec;
    Density rho_tot(false);
    computeDensities(rho_el, rho_tot);
    Timer t_vac;
    mrcpp::ComplexFunction V_vac;
    V_vac.alloc(NUMBER::Real);
    mrcpp::apply(this->apply_prec, V_vac.real(), *poisson, rho_tot.real());
    rho_tot.free(NUMBER::Real);
    print_utils::qmfunction(3, "Vacuum potential", V_vac, t_vac);

    // set up the zero-th iteration potential and gamma, so the first iteration gamma and potentials can be made

    Timer t_gamma;
    if (not this->Vr_n.hasReal()) {
        mrcpp::ComplexFunction gamma_0;
        mrcpp::ComplexFunction V_tot;
        computeGamma(V_vac, gamma_0);
        this->Vr_n = solvePoissonEquation(gamma_0, rho_el);
    }

    // update the potential/gamma before doing anything with them

    Timer t_scrf;
    runMicroIterations(V_vac, rho_el);
    print_utils::qmfunction(3, "Reaction potential", this->Vr_n, t_scrf);
    return this->Vr_n;
}

auto GPESolver::computeEnergies(const Density &rho_el) -> std::tuple<double, double> {
    auto Er_nuc = 0.5 * mrcpp::cplxfunc::dot(this->rho_nuc, this->Vr_n).real();
    auto Er_el = 0.5 * mrcpp::cplxfunc::dot(rho_el, this->Vr_n).real();
    return {Er_el, Er_nuc};
}

void GPESolver::resetComplexFunction(mrcpp::ComplexFunction &function) {
    if (function.hasReal()) function.free(NUMBER::Real);
    if (function.hasImag()) function.free(NUMBER::Imag);
    function.alloc(NUMBER::Real);
}

void GPESolver::printParameters() const {
    nlohmann::json data = {
        {"Method                ", this->solver_name},
        {"Density               ", this->density_type},
        {"Max iterations        ", this->max_iter},
        {"KAIN solver           ", (this->history > 0) ? std::to_string(this->history) : "Off"},
        {"Dynamic threshold     ", (this->dynamic_thrs) ? "On" : "Off"},
    };

    mrcpp::print::separator(3, '~');
    print_utils::json(3, data, false);
    mrcpp::print::separator(3, '~');
}
} // namespace mrchem
