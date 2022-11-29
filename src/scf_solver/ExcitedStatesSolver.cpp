/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "ExcitedStatesSolver.h"
#include "HelmholtzVector.h"
#include "KAIN.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/qmoperator_utils.h"
#include "qmoperators/two_electron/FockBuilder.h"
#include "utils/print_utils.h"

#include <Eigen/Core>
#include <Eigen/StdVector>

using mrcpp::Printer;
using mrcpp::Timer;
using nlohmann::json;

namespace mrchem {

/** @brief Run orbital optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm iterates
 * the Sternheimer response equations in integral form. Common implementation for
 * static and dynamic response. Main points of the algorithm:
 *
 * Pre SCF: setup Helmholtz operators with unperturbed energies
 *
 *  1) Setup perturbed Fock operator
 *  2) For X and Y orbitals do:
 *     a) Apply Helmholtz operator on all orbitals
 *     b) Project out occupied space (1 - rho_0)
 *     c) Compute updates and errors
 *     d) Compute KAIN updates
 *  4) Compute property
 *  5) Check for convergence
 *
 *  Only one state is optimized for now. random guess might help to get more states out
 *
 */
json ExcitedStatesSolver::optimize(Molecule &mol, FockBuilder &F_0, FockBuilder &F_1, int state) {
    Timer t_tot;
    json json_out;
    double err_o = 1.0;
    double err_t = 1.0;

    // Setup KAIN accelerators
    KAIN kain_x(this->history);
    KAIN kain_y(this->history);
    OrbitalVector &Phi_0 = mol.getOrbitals();

    OrbitalVector &X_n = mol.getOrbitalsX(state);

    ComplexMatrix &F_mat_0 = mol.getFockMatrix();

    RankZeroOperator V_0 = F_0.potential();
    bool use_harrison = false; // use harrison´s update scheme for the excitation energies
    bool update_omega = false; // use Kottmann´s update scheme for the excitation energies

    double orb_prec = adjustPrecision(err_o);

    RankZeroOperator V_1 = F_1.potential();
    V_1.setup(orb_prec);

    this->error.push_back(err_t);

    F_0.setup(orb_prec);

    DoubleVector errors_x = DoubleVector::Zero(Phi_0.size());
    DoubleVector errors_y = DoubleVector::Zero(Phi_0.size());

    // orthogonalize all orbitals wrt. the ground state and each other.

    orbital::orthogonalize(this->orth_prec, X_n, Phi_0);
    orbital::orthogonalize(this->orth_prec, X_n);
    // orbital::orthogonalize(this->orth_prec, Y_n, Phi_0);
    // orbital::orthogonalize(this->orth_prec, Y_n);

    auto plevel = Printer::getPrintLevel();

    bool converged = false;
    json_out["cycles"] = {};

    // compute initial omega

    double omega_n = computeOmega(Phi_0, X_n, F_0, V_1, F_mat_0);
    double domega_n = omega_n;

    this->property.push_back(omega_n);

    printParameters(omega_n, F_1.perturbation().name()); // have to change a bit this here
    if (plevel < 1) {
        printConvergenceHeader("Excited state energy");
        if (plevel < 1) printConvergenceRow(0);
    }

    for (auto nIter = 0; (nIter < this->maxIter) or (this->maxIter < 0); nIter++) {
        json json_cycle;
        std::stringstream o_header;
        o_header << "SCF cycle " << nIter;
        mrcpp::print::header(1, o_header.str(), 0, '#');
        mrcpp::print::separator(2, ' ', 1);
        // Initialize SCF cycle
        Timer t_scf, t_lap;
        double orb_prec = adjustPrecision(err_o);
        double helm_prec = getHelmholtzPrec();

        ComplexMatrix F_mat_x = F_mat_0 + omega_n * ComplexMatrix::Identity(Phi_0.size(), Phi_0.size());

        // Setup Helmholtz operators (fixed, based on unperturbed system)
        HelmholtzVector H_x(helm_prec, F_mat_x.real().diagonal());
        auto dot_of_X = orbital::dot(X_n, X_n).sum();

        RankZeroOperator V_1 = F_1.potential();

        if (dynamic and plevel == 1) mrcpp::print::separator(1, '-');

        { // Iterate X orbitals
            // Compute argument: psi_i = sum_i F_0*x_j + (1 - rho_0)V_1(phi_i)
            Timer t_arg;
            mrcpp::print::header(2, "Computing Helmholtz argument");
            t_lap.start();
            V_1.setup(orb_prec);
            OrbitalVector Psi = V_1(Phi_0);
            mrcpp::print::time(2, "Applying V_1", t_lap);
            t_lap.start();
            orbital::orthogonalize(this->orth_prec, Psi, Phi_0);
            mrcpp::print::time(2, "Projecting (1 - rho_0)", t_lap);

            mrcpp::print::footer(2, t_arg, 2);
            if (plevel == 1) mrcpp::print::time(1, "Computing Helmholtz argument", t_arg);

            // Apply Helmholtz operators
            OrbitalVector X_np1 = H_x.apply(V_0, X_n, Psi);
            Psi.clear();
            // Projecting (1 - rho_0)X
            mrcpp::print::header(2, "Projecting occupied space");
            t_lap.start();
            orbital::orthogonalize(this->orth_prec, X_np1, Phi_0);
            orbital::orthogonalize(this->orth_prec, X_np1);

            mrcpp::print::time(2, "Projecting (1 - rho_0)", t_lap);
            mrcpp::print::footer(2, t_lap, 2);

            if (plevel == 1) mrcpp::print::time(1, "Projecting occupied space", t_lap);

            // Compute update and errors
            OrbitalVector dX_n = orbital::add(1.0, X_np1, -1.0, X_n);
            if (update_omega) { domega_n = updateOmega(X_n, X_np1); }
            errors_x = orbital::get_norms(dX_n);

            // Compute KAIN update:
            kain_x.accelerate(orb_prec, X_n, dX_n);

            if (use_harrison) {
                auto V_0_x = V_0(X_n);
                auto left_hand = orbital::add(1.0, V_0_x, 1.0, Psi);
                X_np1 = orbital::add(1.0, X_n, 1.0, dX_n);
                domega_n = -orbital::dot(left_hand, dX_n).sum().real() / orbital::dot(X_np1, X_np1).sum().real();
            }

            Psi.clear();

            // Prepare for next iteration
            X_n = orbital::add(1.0, X_n, 1.0, dX_n);
            // orbital::orthogonalize(this->orth_prec, X_n);

            // Setup perturbed Fock operator (including V_1)
            V_1.clear();
            V_1.setup(orb_prec); // do the x orbitals being uodated change this? it obviously should, but check

            // Compute omega
            mrcpp::print::header(2, "Computing frequency update");
            t_lap.start();

            auto omega_np1 = computeOmega(Phi_0, X_n, F_0, V_1, F_mat_0); /*  computeOmega(Phi_0, X_n, Y_n, F_0, F_1, F_mat_0); */
            domega_n = omega_np1 - omega_n;                               // maybe I should do this before normalization
            omega_n += domega_n;
            mrcpp::print::footer(2, t_lap, 2);
            if (plevel == 1) mrcpp::print::time(1, "Computing frequency update", t_lap);
            this->property.push_back(omega_n);
            X_np1.clear();

            // Save checkpoint file
            if (this->checkpoint) orbital::save_orbitals(X_n, this->chkFileX);
        }

        // Compute errors
        err_o = std::max(errors_x.maxCoeff(), errors_y.maxCoeff());
        err_t = std::sqrt(errors_x.dot(errors_x) + errors_y.dot(errors_y));
        json_cycle["mo_residual"] = err_t;
        // Collect convergence data
        this->error.push_back(err_t);
        double err_w = domega_n;
        converged = checkConvergence(err_o, err_w);
        json_cycle["frequency"] = omega_n;
        json_cycle["frequency_update"] = err_w;

        // Finalize SCF cycle
        if (plevel < 1) printConvergenceRow(nIter + 1);
        printOrbitals(orbital::get_norms(X_n), errors_x, X_n, 1);
        mrcpp::print::separator(1, '-');
        printResidual(err_t, converged);
        mrcpp::print::separator(2, '=', 2);
        printProperty();
        printMemory();
        t_scf.stop();
        json_cycle["wall_time"] = t_scf.elapsed();
        mrcpp::print::footer(1, t_scf, 2, '#');
        mrcpp::print::separator(2, ' ', 2);
        json_out["cycles"].push_back(json_cycle);
        if (converged) {
            mrcpp::print::header(2, "Computing frequency");
            t_lap.start();
            omega_n = computeOmega(Phi_0, X_n, F_0, V_1, F_mat_0); /*  computeOmega(Phi_0, X_n, Y_n, F_0, F_1, F_mat_0); */
            mrcpp::print::footer(2, t_lap, 2);
            if (plevel == 1) mrcpp::print::time(1, "Computing frequency", t_lap);
            this->property.push_back(omega_n);
            V_1.clear();
            break;
        }
        // Clear perturbed Fock operator
        V_1.clear();
    }
    // Compute property
    printConvergence(converged, "Excited states");
    reset();
    json_out["frequency"] = omega_n;
    json_out["wall_time"] = t_tot.elapsed();
    json_out["converged"] = converged;
    return json_out;
}

/** @brief consider only diagonals of the A and S matrices, first for single state. maybe only ok for tda */
double ExcitedStatesSolver::computeOmega(OrbitalVector &Phi, OrbitalVector &X, FockBuilder &F_0, RankZeroOperator &V_1, ComplexMatrix &F_mat_0) {
    // complexvector containing all terms <x_i|x_i>
    auto xi_t_xi_vec = orbital::dot(X, X);
    auto ei_vec = F_mat_0.diagonal(); // could probably append X to Phi_0 for this and then remove them after the computation
    auto sum_ei_xi_t_xi = ei_vec.dot(xi_t_xi_vec).real();
    auto xi_t_F_0_xi_vec = F_0(X, X).trace().real();
    OrbitalVector Psi = V_1(Phi);

    orbital::orthogonalize(this->orth_prec, Psi, Phi);

    auto xi_t_Q_Psi_vec = orbital::dot(X, Psi).sum().real();

    auto F_rr = xi_t_F_0_xi_vec + xi_t_Q_Psi_vec;
    auto X_dot_X = xi_t_xi_vec.sum().real();
    auto omega = (F_rr - sum_ei_xi_t_xi) / X_dot_X;

    return omega;
}

double ExcitedStatesSolver::updateOmega(OrbitalVector &X_n, OrbitalVector &X_np1) {
    auto X_dot_GVX = orbital::dot(X_n, X_np1).sum().real();
    auto GVX_GVX = orbital::dot(X_np1, X_np1).sum().real();
    auto domega_n = X_dot_GVX / GVX_GVX;
    return domega_n;
}

/** @brief Pretty printing of the computed property with update */
void ExcitedStatesSolver::printProperty() const {
    double prop_0(0.0), prop_1(0.0);
    int iter = this->property.size();
    if (iter > 1) prop_0 = this->property[iter - 2];
    if (iter > 0) prop_1 = this->property[iter - 1];

    int w0 = (Printer::getWidth() - 1);
    int w1 = 20;
    int w2 = w0 / 3;
    int w3 = 8;
    int w4 = w0 - w1 - w2 - w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << " ";
    o_head << std::setw(w2) << "Value";
    o_head << std::setw(w4) << "Update";
    o_head << std::setw(w3) << "Done";

    mrcpp::print::separator(2, '=');
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    printUpdate(1, " frequency", prop_1, prop_1 - prop_0, this->propThrs);
    mrcpp::print::separator(2, '=', 2);
}

void ExcitedStatesSolver::printParameters(double omega, const std::string &oper) const {
    std::stringstream o_calc;
    o_calc << "Optimize linear response orbitals";

    std::stringstream o_omega;
    if (this->dynamic) {
        o_omega << std::setprecision(5) << std::fixed << omega << " au";
    } else {
        o_omega << "Static field";
    }

    std::stringstream o_kain;
    if (this->history > 0) {
        o_kain << this->history;
    } else {
        o_kain << "Off";
    }
    std::stringstream o_iter;
    if (this->maxIter > 0) {
        o_iter << this->maxIter;
    } else {
        o_iter << "Off";
    }

    std::stringstream o_prec_0, o_prec_1;
    o_prec_0 << std::setprecision(5) << std::scientific << this->orbPrec[1];
    o_prec_1 << std::setprecision(5) << std::scientific << this->orbPrec[2];

    std::stringstream o_thrs_p;
    if (this->propThrs < 0.0) {
        o_thrs_p << "Off";
    } else {
        o_thrs_p << std::setprecision(5) << std::scientific << this->propThrs;
    }

    std::stringstream o_thrs_o;
    if (this->orbThrs < 0.0) {
        o_thrs_o << "Off";
    } else {
        o_thrs_o << std::setprecision(5) << std::scientific << this->orbThrs;
    }

    std::stringstream o_helm;
    if (this->helmPrec < 0.0) {
        o_helm << "Dynamic";
    } else {
        o_helm << std::setprecision(5) << std::scientific << this->helmPrec;
    }

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation        ", o_calc.str());
    print_utils::text(0, "Frequency          ", o_omega.str());
    print_utils::text(0, "Perturbation       ", oper);
    print_utils::text(0, "Method             ", this->methodName);
    print_utils::text(0, "Relativity         ", this->relativityName);
    print_utils::text(0, "Checkpointing      ", (this->checkpoint) ? "On" : "Off");
    print_utils::text(0, "Max iterations     ", o_iter.str());
    print_utils::text(0, "KAIN solver        ", o_kain.str());
    print_utils::text(0, "Start precision    ", o_prec_0.str());
    print_utils::text(0, "Final precision    ", o_prec_1.str());
    print_utils::text(0, "Helmholtz precision", o_helm.str());
    print_utils::text(0, "Property threshold ", o_thrs_p.str());
    print_utils::text(0, "Orbital threshold  ", o_thrs_o.str());
    mrcpp::print::separator(0, '~', 2);
}

} // namespace mrchem
