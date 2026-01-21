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

#include "GroundStateSolver.h"
#include "HelmholtzVector.h"
#include "KAIN.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/ZoraOperator.h"
#include "qmoperators/two_electron/FockBuilder.h"
#include "qmoperators/two_electron/ReactionOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;
using nlohmann::json;

namespace mrchem {

/** @brief Computes the SCF energy update from last iteration */
double GroundStateSolver::calcPropertyError() const {
    int iter = this->property.size();
    return std::abs(getUpdate(this->property, iter, true));
}

/** @brief Pretty printing of the different contributions to the SCF energy */
void GroundStateSolver::printProperty() const {
    SCFEnergy scf_0, scf_1;
    int iter = this->energy.size();
    if (iter > 1) scf_0 = this->energy[iter - 2];
    if (iter > 0) scf_1 = this->energy[iter - 1];

    double T_0 = scf_0.getKineticEnergy();
    double T_1 = scf_1.getKineticEnergy();
    double V_0 = scf_0.getElectronNuclearEnergy();
    double V_1 = scf_1.getElectronNuclearEnergy();
    double J_0 = scf_0.getElectronElectronEnergy();
    double J_1 = scf_1.getElectronElectronEnergy();
    double K_0 = scf_0.getExchangeEnergy();
    double K_1 = scf_1.getExchangeEnergy();
    double XC_0 = scf_0.getExchangeCorrelationEnergy();
    double XC_1 = scf_1.getExchangeCorrelationEnergy();
    double E_0 = scf_0.getElectronicEnergy();
    double E_1 = scf_1.getElectronicEnergy();
    double N_0 = scf_0.getNuclearEnergy();
    double N_1 = scf_1.getNuclearEnergy();
    double E_eext_0 = scf_0.getElectronExternalEnergy();
    double E_eext_1 = scf_1.getElectronExternalEnergy();
    double E_next_0 = scf_0.getNuclearExternalEnergy();
    double E_next_1 = scf_1.getNuclearExternalEnergy();
    double Er_0 = scf_0.getReactionEnergy();
    double Er_1 = scf_1.getReactionEnergy();
    double Er_el_0 = scf_0.getElectronReactionEnergy();
    double Er_el_1 = scf_1.getElectronReactionEnergy();
    double Er_nuc_0 = scf_0.getNuclearReactionEnergy();
    double Er_nuc_1 = scf_1.getNuclearReactionEnergy();

    bool has_react = (std::abs(Er_el_1) > mrcpp::MachineZero) || (std::abs(Er_nuc_1) > mrcpp::MachineZero);
    bool has_ext = (std::abs(E_eext_1) > mrcpp::MachineZero) || (std::abs(E_next_1) > mrcpp::MachineZero);

    int w0 = (Printer::getWidth() - 1);
    int w1 = 20;
    int w2 = w0 / 3;
    int w3 = 8;
    int w4 = w0 - w1 - w2 - w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << " ";
    o_head << std::setw(w2) << "Energy";
    o_head << std::setw(w4) << "Update";
    o_head << std::setw(w3) << "Done";

    mrcpp::print::separator(2, '=');
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    printUpdate(2, " Kinetic energy  ", T_1, T_1 - T_0, this->propThrs);
    printUpdate(2, " N-E energy      ", V_1, V_1 - V_0, this->propThrs);
    printUpdate(2, " Coulomb energy  ", J_1, J_1 - J_0, this->propThrs);
    printUpdate(2, " Exchange energy ", K_1, K_1 - K_0, this->propThrs);
    printUpdate(2, " X-C energy      ", XC_1, XC_1 - XC_0, this->propThrs);
    if (has_ext) {
        mrcpp::print::separator(2, '-');
        printUpdate(2, " External field (el)  ", E_eext_1, E_eext_1 - E_eext_0, this->propThrs);
        printUpdate(2, " External field (nuc) ", E_next_1, E_next_1 - E_next_0, this->propThrs);
        printUpdate(2, " External field (tot) ", (E_eext_1 + E_next_1), (E_eext_1 + E_next_1) - (E_eext_0 + E_next_0), this->propThrs);
    }
    if (has_react) {
        mrcpp::print::separator(2, '-');
        printUpdate(2, " Reaction energy (el) ", Er_el_1, Er_el_1 - Er_el_0, this->propThrs);
        printUpdate(2, " Reaction energy (nuc) ", Er_nuc_1, Er_nuc_1 - Er_nuc_0, this->propThrs);
        printUpdate(2, " Reaction energy (tot)  ", Er_1, Er_1 - Er_0, this->propThrs);
    }
    mrcpp::print::separator(2, '-');
    printUpdate(2, " Electronic energy", E_1, E_1 - E_0, this->propThrs);
    printUpdate(2, " Nuclear energy   ", N_1, N_1 - N_0, this->propThrs);
    mrcpp::print::separator(2, '-');
    printUpdate(1, " Total energy     ", E_1 + N_1, (E_1 + N_1) - (E_0 + N_0), this->propThrs);
    mrcpp::print::separator(2, '=', 2);
}

void GroundStateSolver::printParameters(const std::string &calculation) const {
    std::stringstream o_iter;
    if (this->maxIter > 0) {
        o_iter << this->maxIter;
    } else {
        o_iter << "Off";
    }

    std::stringstream o_kain;
    if (this->history > 0) {
        o_kain << this->history;
    } else {
        o_kain << "Off";
    }

    std::stringstream o_loc;
    if (this->localize) {
        if (this->rotation == 0) {
            o_loc << "First two iterations";
        } else if (this->rotation == 1) {
            o_loc << "Every iteration";
        } else {
            o_loc << "Every " << this->rotation << " iterations";
        }
    } else {
        o_loc << "Off";
    }

    std::stringstream o_diag;
    if (not this->localize) {
        if (this->rotation == 0) {
            o_diag << "First two iterations";
        } else if (this->rotation == 1) {
            o_diag << "Every iteration";
        } else {
            o_diag << "Every " << this->rotation << " iterations";
        }
    } else {
        o_diag << "Off";
    }

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

    std::stringstream o_prec_0, o_prec_1;
    o_prec_0 << std::setprecision(5) << std::scientific << this->orbPrec[1];
    o_prec_1 << std::setprecision(5) << std::scientific << this->orbPrec[2];

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation        ", calculation);
    print_utils::text(0, "Method             ", this->methodName);
    print_utils::text(0, "Relativity         ", this->relativityName);
    print_utils::text(0, "Environment        ", this->environmentName);
    print_utils::text(0, "External fields    ", this->externalFieldName);
    print_utils::text(0, "Checkpointing      ", (this->checkpoint) ? "On" : "Off");
    print_utils::text(0, "Max iterations     ", o_iter.str());
    print_utils::text(0, "KAIN solver        ", o_kain.str());
    print_utils::text(0, "Localization       ", o_loc.str());
    print_utils::text(0, "Diagonalization    ", o_diag.str());
    print_utils::text(0, "Start precision    ", o_prec_0.str());
    print_utils::text(0, "Final precision    ", o_prec_1.str());
    print_utils::text(0, "Helmholtz precision", o_helm.str());
    print_utils::text(0, "Energy threshold   ", o_thrs_p.str());
    print_utils::text(0, "Orbital threshold  ", o_thrs_o.str());
    mrcpp::print::separator(0, '~', 2);
}

/** @brief Reset accumulated data */
void GroundStateSolver::reset() {
    SCFSolver::reset();
    this->energy.clear();
}

/** @brief Run orbital optimization
 *
 * @param mol: Molecule to optimize
 * @param F: FockBuilder defining the SCF equations
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm computes
 * the Fock matrix explicitly using the kinetic energy operator, and uses a KAIN
 * accelerator to improve convergence. Diagonalization or localization may be performed
 * during the SCF iterations. Main points of the algorithm:
 *
 * Pre SCF: setup Fock operator and compute Fock matrix
 *
 *  1) Diagonalize/localize orbitals
 *  2) Compute current SCF energy
 *  3) Apply Helmholtz operator on all orbitals
 *  4) Orthonormalize orbitals (Löwdin)
 *  5) Compute orbital updates
 *  6) Compute KAIN update
 *  7) Compute errors and check for convergence
 *  8) Add orbital updates
 *  9) Orthonormalize orbitals (Löwdin)
 * 10) Setup Fock operator
 * 11) Compute Fock matrix
 *
 */
json GroundStateSolver::optimize(Molecule &mol, FockBuilder &F) {
    printParameters("Optimize ground state orbitals");

    Timer t_tot;
    json json_out;

    SCFEnergy &E_n = mol.getSCFEnergy();
    const Nuclei &nucs = mol.getNuclei();
    OrbitalVector &Phi_n = mol.getOrbitals();
    ComplexMatrix &F_mat = mol.getFockMatrix();

    auto scaling = std::vector<double>(Phi_n.size(), 1.0);
    KAIN kain(this->history, 0, false, scaling);

    DoubleVector errors = DoubleVector::Ones(Phi_n.size());
    double err_o = errors.maxCoeff();
    double err_t = errors.norm();

    this->error.push_back(err_t);
    this->energy.push_back(E_n);
    this->property.push_back(E_n.getTotalEnergy());

    auto plevel = Printer::getPrintLevel();
    if (plevel < 1) {
        printConvergenceHeader("Total energy");
        printConvergenceRow(0);
    }

    // Initialize Resolvent (Attention: HelmholtzVector = -2 Helmholtz)
    HelmholtzVector Minus_2_Resolvent(getHelmholtzPrec(), Eigen::VectorXd::Constant(Phi_n.size(), -0.5));
    ResolventVector         Resolvent(getHelmholtzPrec(), Eigen::VectorXd::Constant(Phi_n.size(), -1.0));

    // Parameters for line search
    double alpha = 1.0;                          // current step size (adaptive across iterations)
    const double beta = 0.5;                     // shrink factor (0 < beta < 1)
    const double gamma = 1.4;                    // growth factor (>1)
    const double armijo_parameter = 1e-4;        // Armijo parameter
    const double alpha_min = 1e-12;              // safeguard lower bound
    const double alpha_max = 10.0;               // safeguard upper bound

    // Parameters for restarting and momentum
    int last_restart_iter = -1; //0
    const int restart_cooldown = 4;     // no restarts within these many iterations of previous restart
    const double eta_powell = 0.3;         // Powell threshold (tune 0.1..0.3)
    const double polak_max = 5.0;          // cap on beta (safeguard)

    OrbitalVector direction;
    OrbitalVector previous_grad_E;
    OrbitalVector previous_preconditioned_grad_E;

    double previous_h1_inner_product_preconditioned_grad_E_grad_E;

    int nIter = 0;
    bool converged = false;
    json_out["cycles"] = {};
    while (nIter++ < this->maxIter or this->maxIter < 0) {
        json json_cycle;
        std::stringstream o_header;
        o_header << "SCF cycle " << nIter;
        mrcpp::print::header(1, o_header.str(), 0, '#');
        mrcpp::print::separator(2, ' ', 1);

        // Initialize SCF cycle
        Timer t_scf;
        double orb_prec = adjustPrecision(err_o);
        double helm_prec = getHelmholtzPrec();
        if (nIter < 2) {
            if (F.getReactionOperator() != nullptr) F.getReactionOperator()->updateMOResidual(err_t);
            F.setup(orb_prec);
        }

        // Init Helmholtz operator
        HelmholtzVector H(helm_prec, F_mat.real().diagonal());
        ComplexMatrix L_mat = H.getLambdaMatrix();

        //std::cout << "Fock matrix:" << std::endl << F_mat.real() << std::endl;
        //std::cout << "Lambda matrix:" << std::endl << L_mat.real() << std::endl;
        
        // Apply Helmholtz operator
        OrbitalVector Psi = F.buildHelmholtzArgument(orb_prec, Phi_n, F_mat, L_mat);
        OrbitalVector Phi_np1 = H(Psi);

        Psi = F.potential()(Phi_n);
        OrbitalVector grad_E = orbital::add(1.0, Phi_n, -2.0, Psi);
        OrbitalVector grad_E1= orbital::add(1.0, Psi, -0.5, Phi_n);
        Psi.clear();
        grad_E = Minus_2_Resolvent(grad_E);
        grad_E1 = Resolvent(grad_E1);
        grad_E = orbital::add(2.0, Phi_n, 1.0, grad_E);
        grad_E1 = orbital::add(2.0, Phi_n, 4.0, grad_E1);

        OrbitalVector Minus_2_Resolvent_Phi = Minus_2_Resolvent(Phi_n);
        ComplexMatrix B_proj = -0.5 * orbital::calc_overlap_matrix(Minus_2_Resolvent_Phi, Phi_n);
        ComplexMatrix C_proj_complex = orbital::calc_overlap_matrix(grad_E, Phi_n);
        DoubleMatrix C_proj_sym = (C_proj_complex.real() + C_proj_complex.real().transpose()) * 0.5;
        DoubleMatrix B_proj_real = (B_proj.real() + B_proj.real().transpose()) * 0.5;
        DoubleMatrix A_proj = mrchem::math_utils::solve_symmetric_sylvester(B_proj_real, 2.0 * C_proj_sym);
        DoubleMatrix minus_half_A_proj = -0.5 * A_proj;
        
        
        OrbitalVector AR_Phi = orbital::rotate(Minus_2_Resolvent_Phi, minus_half_A_proj);
        grad_E = orbital::add(1.0, grad_E, -1.0, AR_Phi);
        //auto grad_E_norm = orbital::get_norms(grad_E).maxCoeff();
        //std::cout << "norm(grad_E) = " << grad_E_norm << std::endl;
        //std::cout << "--------------------------------------" << std::endl;

        OrbitalVector Resolvent_Phi = Resolvent(Phi_n);
        ComplexMatrix B_proj1 = orbital::calc_overlap_matrix(Resolvent_Phi, Phi_n);
        ComplexMatrix C_proj_complex1 = orbital::calc_overlap_matrix(grad_E1, Phi_n);
        DoubleMatrix C_proj_sym1 = C_proj_complex1.real() + C_proj_complex1.real().transpose();
        DoubleMatrix B_proj_real1 = (B_proj1.real() + B_proj1.real().transpose()) * 0.5;
        DoubleMatrix A_proj1 = mrchem::math_utils::solve_symmetric_sylvester(B_proj_real1, C_proj_sym1);

        OrbitalVector AR_Phi1 = orbital::rotate(Resolvent_Phi, A_proj1);
        grad_E1 = orbital::add(1.0, grad_E1, -1.0, AR_Phi1);
        auto grad_E1_norm = orbital::get_norms(grad_E1).maxCoeff();
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "norm(grad_E1) = " << grad_E1_norm << std::endl;
        std::cout << "--------------------------------------" << std::endl;

        //auto grad_E_error = orbital::add(1.0, grad_E, -1.0, grad_E1);
        //auto grad_E_error_norm = orbital::get_norms(grad_E_error).maxCoeff();
        //std::cout << "--------------------------------------" << std::endl;
        //std::cout << "norm(grad_E - grad_E1) = " << grad_E_error_norm << std::endl;
        //std::cout << "--------------------------------------" << std::endl;


        

        F.clear();

        // ==============================
        // Preconditioning
        // ==============================

        A_proj = A_proj1;

        OrbitalVector preconditioned_grad_E = grad_E;

        // Diagonalize A_proj
        Eigen::SelfAdjointEigenSolver<DoubleMatrix> eigensolver(A_proj);
        if (eigensolver.info() != Eigen::Success) {
            MSG_ABORT("Eigen-decomposition of A_proj failed");
        }

        Eigen::VectorXd sigma_A_proj = eigensolver.eigenvalues();
        DoubleMatrix U_A_proj = eigensolver.eigenvectors();

        // Check norm(A - 4F) tends to zero
        std::cout << "norm(A_proj - 4 * F_mat.real()) = " 
                  << (A_proj - 4.0 * F_mat.real()).norm() << std::endl;

        
        
        if (sigma_A_proj.maxCoeff() < 0.0) {

            Eigen::VectorXd orbital_energy = 0.5 * sigma_A_proj;
            Eigen::MatrixXd one_plus_orbital_energy = (Eigen::VectorXd::Ones(orbital_energy.size()) + orbital_energy).asDiagonal();

            ResolventVector Resolvent_mu( getHelmholtzPrec(), orbital_energy );

            preconditioned_grad_E = orbital::rotate(preconditioned_grad_E, U_A_proj.transpose());

            auto temp = Resolvent_mu(preconditioned_grad_E);
            temp = orbital::rotate(temp, one_plus_orbital_energy);
            preconditioned_grad_E = orbital::add( 0.5, preconditioned_grad_E, 0.5, temp );

            preconditioned_grad_E = orbital::rotate(preconditioned_grad_E, U_A_proj);
        }

        // Check norm of preconditioned gradient
        auto &nabla = F.momentum();
        nabla.setup(orb_prec);

        auto preconditioned_grad_E_norm = orbital::h1_norm(preconditioned_grad_E, nabla);
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "norm(preconditioned_grad_E) = " << preconditioned_grad_E_norm << std::endl;
        std::cout << "--------------------------------------" << std::endl;

        C_proj_complex1 = orbital::calc_overlap_matrix(preconditioned_grad_E, Phi_n);
        C_proj_sym1 = C_proj_complex1.real() + C_proj_complex1.real().transpose();
        A_proj1 = mrchem::math_utils::solve_symmetric_sylvester(B_proj_real1, C_proj_sym1);
        AR_Phi1 = orbital::rotate(Resolvent_Phi, A_proj1);
        preconditioned_grad_E = orbital::add(1.0, preconditioned_grad_E, -1.0, AR_Phi1);

        grad_E1_norm = orbital::h1_norm(preconditioned_grad_E, nabla);
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "norm(preconditioned_grad_E and projected) = " << grad_E1_norm << std::endl;
        std::cout << "--------------------------------------" << std::endl;
        
        //# Necessary for Grassmann: 
        //preconditioned_grad_E = project_to_horizontal(preconditioned_grad_E, Phi)
        
        // ==============================
        // End Preconditioning

        auto grad_E_norm = orbital::h1_norm(grad_E, nabla);
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "norm(grad_E) = " << grad_E_norm << std::endl;
        std::cout << "--------------------------------------" << std::endl;

        // Safeguard: if not descent direction, skip preconditioning
        double h1_inner_product_preconditioned_grad_E_grad_E = orbital::h1_inner_product(preconditioned_grad_E, grad_E, nabla);
        if (h1_inner_product_preconditioned_grad_E_grad_E <= 0.0) {
            std::cout << "Preconditioning skipped (not a descent direction): " << h1_inner_product_preconditioned_grad_E_grad_E << std::endl;
            preconditioned_grad_E = grad_E;
            h1_inner_product_preconditioned_grad_E_grad_E = grad_E_norm * grad_E_norm;
        }

        // ======================================================
        // Conjugate-gradient direction (H1, Polak-Ribière)
        // ======================================================
        
        if (nIter == 1) {
            // First iteration: steepest descent
            direction = orbital::add(-1.0, preconditioned_grad_E, 0.0, preconditioned_grad_E);
        }
        else {
            // Polak–Ribière coefficient
            OrbitalVector diff_pc_grad = orbital::add(1.0, preconditioned_grad_E, -1.0, previous_preconditioned_grad_E);
            double polak_ribiere = orbital::h1_inner_product(diff_pc_grad, grad_E, nabla);
            polak_ribiere = polak_ribiere / (previous_h1_inner_product_preconditioned_grad_E_grad_E + mrcpp::MachineZero);
            polak_ribiere = std::max(0.0, polak_ribiere);
            polak_ribiere = std::min(polak_ribiere, polak_max);
            std::cout << "Polak-Ribière coefficient = " << polak_ribiere << std::endl;

            // Project previous direction to tangent space
            ComplexMatrix C_proj_dir = orbital::calc_overlap_matrix(direction, Phi_n);
            DoubleMatrix C_proj_dir_sym = (C_proj_dir.real() + C_proj_dir.real().transpose()) * 0.5;
            DoubleMatrix A_proj_dir = mrchem::math_utils::solve_symmetric_sylvester(B_proj_real1, C_proj_dir_sym);

            OrbitalVector projected_direction = orbital::rotate(Resolvent_Phi, A_proj_dir);
            projected_direction = orbital::add(1.0, direction, -1.0, projected_direction);
            // Grassmann horizontal projection (optional but recommended)
            //projected_direction = orbital::project_to_horizontal(projected_direction, Phi_n);

            direction = orbital::add(polak_ribiere, projected_direction, -1.0, preconditioned_grad_E);

            // ---------- Robust restart checks ----------
            bool do_restart = false;
            auto reason = "no reason";

            // (1) Descent check
            double descent = orbital::h1_inner_product(direction, grad_E, nabla);

            if (descent >= 0.0) {
                do_restart = true;
                reason = "non-descent";
            }

            // (2) Powell restart
            else if ((nIter - last_restart_iter) > restart_cooldown) {
                double inner =
                    orbital::h1_inner_product(grad_E, previous_preconditioned_grad_E, nabla);

                double ref =
                    orbital::h1_inner_product(previous_preconditioned_grad_E, previous_grad_E, nabla);

                if (inner >= eta_powell * ref)
                    do_restart = true;
                    reason = "powell";
            }

            if (do_restart) {
                std::cout << "Powell/guarded restart at iteration_index " << nIter << " (reason: " << reason << ")" << std::endl;
                direction = orbital::add(-1.0, preconditioned_grad_E, 0.0, preconditioned_grad_E);
                last_restart_iter = nIter;
            }
            // ------------------------------------------
        }



        previous_preconditioned_grad_E = preconditioned_grad_E;
        previous_grad_E = grad_E;
        previous_h1_inner_product_preconditioned_grad_E_grad_E = h1_inner_product_preconditioned_grad_E_grad_E;

        // Backtracking line search
        OrbitalVector Phi_backup = orbital::deep_copy(Phi_n);
        auto alpha_trial = alpha;
        int count = 0;
        while (true) {
            count += 1;
            // Retraction to Stiefel is Lowdin based:
            Phi_n = orbital::add(1.0, Phi_n, alpha_trial, direction);
            // Orthonormalization updates F_mat as a side effect?!
            orbital::orthonormalize(orb_prec, Phi_n, F_mat);
            // Compute Fock matrix and energy
            F.setup(orb_prec);
            //F_mat = F(Phi_candidate, Phi_candidate);
            E_n = F.trace(Phi_n, nucs);
            std::cout << "Candidate Energy: " << E_n.getTotalEnergy() << std::endl;

            break;
        }
        F.clear();

        // Orthonormalize
        orbital::orthonormalize(orb_prec, Phi_np1, F_mat);

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);
        Phi_np1.clear();

        kain.accelerate(orb_prec, Phi_n, dPhi_n);

        // Compute errors
        errors = orbital::get_norms(dPhi_n);
        err_o = errors.maxCoeff();
        err_t = errors.norm();
        json_cycle["mo_residual"] = err_t;

        // Update orbitals
        Phi_n = orbital::add(1.0, Phi_n, 1.0, dPhi_n);
        dPhi_n.clear();

        orbital::orthonormalize(orb_prec, Phi_n, F_mat);

        // Compute Fock matrix and energy
        if (F.getReactionOperator() != nullptr) F.getReactionOperator()->updateMOResidual(err_t);
        F.setup(orb_prec);
        F_mat = F(Phi_n, Phi_n);
        E_n = F.trace(Phi_n, nucs);

        //h1 = orbital::h1_inner_product(Phi_n, Phi_n, nabla);
        //std::cout << "Should be the same: " << h1 << " = " << orbital::h1_norm(Phi_n, nabla) * orbital::h1_norm(Phi_n, nabla) << std::endl;
        //h1 = h1 - orbital::l2_inner_product(Phi_n, Phi_n);
        //h1 = h1 * 0.5;
        //std::cout << "Kinetic energy through h1_inner_product() = " << h1 << std::endl;
// ==============================================================
        // Collect convergence data
        this->error.push_back(err_t);
        this->energy.push_back(E_n);
        this->property.push_back(E_n.getTotalEnergy());
        auto err_p = calcPropertyError();
        converged = checkConvergence(err_o, err_p);

        json_cycle["energy_terms"] = E_n.json();
        json_cycle["energy_total"] = E_n.getTotalEnergy();
        json_cycle["energy_update"] = err_p;

        // Rotate orbitals
        if (needLocalization(nIter, converged)) {
            ComplexMatrix U_mat = orbital::localize(orb_prec, Phi_n, F_mat);
            F.rotate(U_mat);
            kain.clear();
        } else if (needDiagonalization(nIter, converged)) {
            ComplexMatrix U_mat = orbital::diagonalize(orb_prec, Phi_n, F_mat);
            F.rotate(U_mat);
            kain.clear();
        }

        // Save checkpoint file
        if (this->checkpoint) orbital::save_orbitals(Phi_n, this->chkFile);

        // Finalize SCF cycle
        if (plevel < 1) printConvergenceRow(nIter);
        printOrbitals(F_mat.real().diagonal(), errors, Phi_n, 0);
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
        if (converged) break;
    }

    F.clear();
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);

    printConvergence(converged, "Total energy");
    reset();

    Timer t_eps;
    mrcpp::print::header(1, "Computing orbital energies");
    OrbitalEnergies &eps = mol.getOrbitalEnergies();
    eps.getOccupation() = orbital::get_occupations(Phi_n);
    eps.getEpsilon() = orbital::calc_eigenvalues(Phi_n, F_mat);
    eps.getSpin() = orbital::get_spins(Phi_n);
    mrcpp::print::footer(1, t_eps, 2);

    json_out["wall_time"] = t_tot.elapsed();
    json_out["converged"] = converged;
    return json_out;
}

/** @brief Test if orbitals needs localization
 *
 * @param nIter: current iteration number
 *
 * This check is based on the "localize" and "rotation" parameters, where the latter
 * tells how oftern (in terms of iterations) the orbitals should be rotated.
 */
bool GroundStateSolver::needLocalization(int nIter, bool converged) const {
    bool loc = false;
    if (not this->localize) {
        loc = false;
    } else if (nIter <= 2 or converged) {
        loc = true;
    } else if (this->rotation == 0) {
        loc = false;
    } else if (nIter % this->rotation == 0) {
        loc = true;
    }
    return loc;
}

/** @brief Test if orbitals needs diagonalization
 *
 * @param nIter: current iteration number
 *
 * This check is based on the "localize" and "rotation" parameters, where the latter
 * tells how oftern (in terms of iterations) the orbitals should be rotated.
 */
bool GroundStateSolver::needDiagonalization(int nIter, bool converged) const {
    bool diag = false;
    if (this->localize) {
        diag = false;
    } else if (nIter <= 2 or converged) {
        diag = true;
    } else if (this->rotation == 0) {
        diag = false;
    } else if (nIter % this->rotation == 0) {
        diag = true;
    }
    return diag;
}

} // namespace mrchem
