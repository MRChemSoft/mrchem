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
        o_kain << "Grassmann";
    } else {
        o_kain << "Stiefel";
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

/** @brief Run conjugate orbital optimization
 *
 * @param mol: Molecule to optimize
 * @param F: FockBuilder defining the SCF equations
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
    ResolventVector Resolvent(getHelmholtzPrec(), Eigen::VectorXd::Constant(Phi_n.size(), -1.0));

    // Parameters for line search
    double alpha = 1.0;                          // current step size (adaptive across iterations)
    const double beta = 0.5;                     // shrink factor (0 < beta < 1)
    const double gamma = 1.4;                    // growth factor (>1)
    const double armijo_parameter = 1e-4;        // Armijo parameter
    const double alpha_min = 1e-10;              // safeguard lower bound
    const double alpha_max = 10.0;               // safeguard upper bound

    // Parameters for restarting and momentum
    int last_restart_iter = 0;
    const int restart_cooldown = 4;        // no restarts within these many iterations of previous restart
    const double eta_powell = 0.3;         // Powell threshold (tune 0.1..0.3)
    const double polak_max = 5.0;          // cap on beta (safeguard)

    OrbitalVector direction = orbital::param_copy(Phi_n);
    OrbitalVector previous_grad_E = orbital::param_copy(Phi_n);
    OrbitalVector previous_preconditioned_grad_E = orbital::param_copy(Phi_n);

    double previous_h1_inner_product_preconditioned_grad_E_grad_E = 0.0;
    int rejectness_count = 0;
    std::vector<double> grad_E_array;
        
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
        if (nIter < 2) {
            if (F.getReactionOperator() != nullptr) F.getReactionOperator()->updateMOResidual(err_t);
            F.setup(orb_prec);
        }

        // Calculate Euclidian gradient starting with the potential evaluation
        OrbitalVector grad_E = orbital::param_copy(Phi_n);
        grad_E = F.potential()(Phi_n);
        grad_E.distribute();

        for (auto &phi_i : grad_E)
        {
            if (mrcpp::mpi::my_func(phi_i))
                phi_i.crop(orb_prec);
        }
        //MPI_Barrier(mrcpp::mpi::comm_wrk);
        mrcpp::print::separator(0, '-');


        // Compute Laplacian of Phi_n and the quantity (1 - Laplacian)grad_E.
        auto &nabla = F.momentum();
        OrbitalVector dx_Phi = nabla[0](Phi_n);
        OrbitalVector dy_Phi = nabla[1](Phi_n);
        OrbitalVector dz_Phi = nabla[2](Phi_n);
        dx_Phi = nabla[0](dx_Phi);
        dy_Phi = nabla[1](dy_Phi);
        dz_Phi = nabla[2](dz_Phi);
        OrbitalVector one_minus_laplacian_grad_E = orbital::param_copy(Phi_n);
        one_minus_laplacian_grad_E = orbital::add(4.0, grad_E, -2.0, dx_Phi);
        one_minus_laplacian_grad_E = orbital::add(1.0, one_minus_laplacian_grad_E, -2.0, dy_Phi);
        one_minus_laplacian_grad_E = orbital::add(1.0, one_minus_laplacian_grad_E, -2.0, dz_Phi);
        
        OrbitalVector one_minus_laplacian_Phi = orbital::param_copy(Phi_n);
        one_minus_laplacian_Phi = orbital::add(1.0, Phi_n, -1.0, dx_Phi);
        one_minus_laplacian_Phi = orbital::add(1.0, one_minus_laplacian_Phi, -1.0, dy_Phi);
        one_minus_laplacian_Phi = orbital::add(1.0, one_minus_laplacian_Phi, -1.0, dz_Phi);
        
        for (auto &phi_i : one_minus_laplacian_Phi)
        {
            if (mrcpp::mpi::my_func(phi_i))
                phi_i.crop(orb_prec);
        }

        dx_Phi.clear();
        dy_Phi.clear();
        dz_Phi.clear();

        // Evaluate resolvent and its quadratic form
        OrbitalVector Resolvent_Phi = Resolvent(Phi_n);
        ComplexMatrix B_proj = orbital::calc_overlap_matrix(Resolvent_Phi, Phi_n);
        DoubleMatrix B_proj_real = (B_proj.real() + B_proj.real().transpose()) * 0.5;

        ComplexMatrix C_proj_complex = orbital::calc_overlap_matrix(one_minus_laplacian_grad_E, Resolvent_Phi);
        DoubleMatrix C_proj_sym = C_proj_complex.real() + C_proj_complex.real().transpose();
        DoubleMatrix A_proj = mrchem::math_utils::solve_symmetric_sylvester(B_proj_real, C_proj_sym);


        // Check norm(A - 4F) tends to zero
        println(0, "norm(A_proj - 4F) = " << (A_proj - 4.0 * F_mat.real()).norm());
        

        OrbitalVector AR_Phi = orbital::rotate(Phi_n, A_proj);
        one_minus_laplacian_grad_E = orbital::add(1.0, one_minus_laplacian_grad_E, -1.0, AR_Phi);
        AR_Phi.clear();

        for (auto &phi_i : one_minus_laplacian_grad_E)
        {
            if (mrcpp::mpi::my_func(phi_i))
                phi_i.crop(orb_prec);
        }

        // Alternative gradient evaluation
        grad_E = Resolvent(one_minus_laplacian_grad_E);
        
        for (auto &phi_i : grad_E)
        {
            if (mrcpp::mpi::my_func(phi_i))
                phi_i.crop(orb_prec);
        }

        // Check norm of gradient: grad_E = grad_E1 up to numerical noise
        auto grad_E_norm = orbital::get_norms(grad_E).norm();
        println(0, "L2-n(grad_E) = " << grad_E_norm);
        
        grad_E_norm = orbital::l2_inner_product(grad_E, one_minus_laplacian_grad_E);
        grad_E_norm = std::sqrt(std::abs(grad_E_norm));
        println(0, "norm(grad_E) = " << grad_E_norm);
        grad_E_array.push_back(grad_E_norm);
        
        
        // ==============================
        // Preconditioning
        // ==============================

        OrbitalVector preconditioned_grad_E = grad_E;

        // Diagonalize A_proj
        Eigen::SelfAdjointEigenSolver<DoubleMatrix> eigensolver(A_proj);
        if (eigensolver.info() != Eigen::Success) {
            MSG_ABORT("Eigen-decomposition of A_proj failed");
        }

        Eigen::VectorXd sigma_A_proj = eigensolver.eigenvalues();
        DoubleMatrix U_A_proj = eigensolver.eigenvectors();

        if (sigma_A_proj.maxCoeff() >= 0.0)
            MSG_ABORT("Non-negative eigenvalue in A_proj, preconditioning not possible");
        auto sigma_A_proj_inv = - sigma_A_proj.cwiseInverse();
        auto sigma0 = sigma_A_proj_inv.minCoeff();
        sigma0 = std::min(sigma0, 0.5);
        auto sigma1 = sigma_A_proj_inv.maxCoeff();
        sigma1 = std::max(sigma1, 0.5);
        println(0, "Sigma0 = " << sigma0);
        println(0, "Sigma1 = " << sigma1);
        
        auto lower_preconditioning_boundary = sigma0 * grad_E_norm * grad_E_norm;
        auto upper_preconditioning_boundary = sigma1 * grad_E_norm * grad_E_norm;
        
        
        Eigen::VectorXd orbital_energy = 0.5 * sigma_A_proj;
        Eigen::MatrixXd one_plus_orbital_energy = (Eigen::VectorXd::Ones(orbital_energy.size()) + orbital_energy).asDiagonal();
        ResolventVector Resolvent_mu( getHelmholtzPrec(), orbital_energy );
        preconditioned_grad_E = orbital::rotate(preconditioned_grad_E, U_A_proj.transpose());
        auto temp = Resolvent_mu(preconditioned_grad_E);
        temp = orbital::rotate(temp, one_plus_orbital_energy);
        preconditioned_grad_E = orbital::add( 0.5, preconditioned_grad_E, 0.5, temp );
        //temp.clear();
        preconditioned_grad_E = orbital::rotate(preconditioned_grad_E, U_A_proj);
        
        for (auto &phi_i : preconditioned_grad_E)
        {
            if (mrcpp::mpi::my_func(phi_i))
            {
                phi_i.crop(orb_prec);
            }
        }

        
        C_proj_complex = orbital::calc_overlap_matrix(preconditioned_grad_E, Phi_n);
        C_proj_sym = C_proj_complex.real() + C_proj_complex.real().transpose();
        A_proj = mrchem::math_utils::solve_symmetric_sylvester(B_proj_real, C_proj_sym);
        AR_Phi = orbital::rotate(Resolvent_Phi, A_proj);
        preconditioned_grad_E = orbital::add(1.0, preconditioned_grad_E, -1.0, AR_Phi);
        AR_Phi.clear();

        
        // Necessary for Grassmann: 
        if (this->history > 0)
        {
            preconditioned_grad_E = orbital::project_to_horizontal(preconditioned_grad_E, Phi_n, nabla);
        }


        double h1_inner_product_preconditioned_grad_E_grad_E = orbital::l2_inner_product(preconditioned_grad_E, one_minus_laplacian_grad_E);
        println(0, "lower_preconditioning_boundary = " << lower_preconditioning_boundary);
        println(0, "product(preconditioned_grad_E, grad_E, 1) = " << h1_inner_product_preconditioned_grad_E_grad_E);
        println(0, "upper_preconditioning_boundary = " << upper_preconditioning_boundary);
        if
        (
            lower_preconditioning_boundary > h1_inner_product_preconditioned_grad_E_grad_E
            ||
            h1_inner_product_preconditioned_grad_E_grad_E > upper_preconditioning_boundary
        )
        {
            println(0, "Preconditioning skipped (noise dominates gradient)");
            preconditioned_grad_E = grad_E;
            h1_inner_product_preconditioned_grad_E_grad_E = grad_E_norm * grad_E_norm;
        }
        // End Preconditioning
        // ==============================
        
        
        // ======================================================
        // Conjugate-gradient direction (H1, Polak-Ribière)
        // ======================================================
        
        double descent_directional_product;

        //if (true) {
        if (nIter == 1) {
            // First iteration: steepest descent
            direction = orbital::add(-1.0, preconditioned_grad_E, 0.0, preconditioned_grad_E);
            descent_directional_product = - h1_inner_product_preconditioned_grad_E_grad_E;
        }
        else {
            // Polak–Ribière coefficient
            OrbitalVector diff_pc_grad = orbital::add(1.0, preconditioned_grad_E, -1.0, previous_preconditioned_grad_E);
            double polak_ribiere = orbital::l2_inner_product(diff_pc_grad, one_minus_laplacian_grad_E);
            polak_ribiere = polak_ribiere / (previous_h1_inner_product_preconditioned_grad_E_grad_E + mrcpp::MachineZero);
            polak_ribiere = std::max(0.0, polak_ribiere);
            polak_ribiere = std::min(polak_ribiere, polak_max);
            println(0, "Polak-Ribiere coefficient = " << polak_ribiere);
            bool polak_ribiere_is_small = (polak_ribiere <= mrcpp::MachineZero);

            if (not polak_ribiere_is_small)
            {
                // Project previous direction to tangent space
                ComplexMatrix C_proj_dir = orbital::calc_overlap_matrix(direction, Phi_n);
                DoubleMatrix C_proj_dir_sym = (C_proj_dir.real() + C_proj_dir.real().transpose()) * 0.5;
                DoubleMatrix A_proj_dir = mrchem::math_utils::solve_symmetric_sylvester(B_proj_real, C_proj_dir_sym);

                OrbitalVector projected_direction = orbital::rotate(Resolvent_Phi, A_proj_dir);
                projected_direction = orbital::add(1.0, direction, -1.0, projected_direction);
                // Necessary for Grassmann:
                if (this->history > 0)
                    projected_direction = orbital::project_to_horizontal(projected_direction, Phi_n, nabla);

                direction = orbital::add(polak_ribiere, projected_direction, -1.0, preconditioned_grad_E);
            }
            else
            {
                direction = orbital::add(-1.0, preconditioned_grad_E, 0.0, preconditioned_grad_E);
                descent_directional_product = - h1_inner_product_preconditioned_grad_E_grad_E;
            }


            // ---------- Robust restart checks ----------
            bool do_restart = false;
            auto reason = "no reason";

            // (1) Descent check
            double descent;
            //if (not polak_ribiere_is_small) descent = orbital::h1_inner_product(direction, grad_E, nabla);
            if (not polak_ribiere_is_small) descent = orbital::l2_inner_product(direction, one_minus_laplacian_grad_E);
            else descent = descent_directional_product;
            
            if (descent >= 0.0) {
                do_restart = true;
                reason = "non-descent";
            }

            // (2) Powell restart
            else if ((nIter - last_restart_iter) > restart_cooldown) {
                double inner =
                    //orbital::h1_inner_product(grad_E, previous_preconditioned_grad_E, nabla);
                    orbital::l2_inner_product(one_minus_laplacian_grad_E, previous_preconditioned_grad_E);

                if (inner >= eta_powell * previous_h1_inner_product_preconditioned_grad_E_grad_E)
                {
                    do_restart = true;
                    reason = "powell";
                }
            }

            if (do_restart) {
                println(0, "Powell/guarded restart (reason: " << reason << ")");
                direction = orbital::add(-1.0, preconditioned_grad_E, 0.0, preconditioned_grad_E);
                descent_directional_product = - h1_inner_product_preconditioned_grad_E_grad_E;
                last_restart_iter = nIter;
            }
            else descent_directional_product = descent;
            // ------------------------------------------
        }



        previous_preconditioned_grad_E = preconditioned_grad_E;
        previous_grad_E = grad_E;
        previous_h1_inner_product_preconditioned_grad_E_grad_E = h1_inner_product_preconditioned_grad_E_grad_E;

        // WARNING: FockBuilder implicitly depends on mol.getOrbitals()
        // Orbital swapping must be scoped and restored
        OrbitalVector Phi_backup = orbital::deep_copy(Phi_n);
        auto Energy = this->property.back();

        // Backtracking line search
        auto alpha_trial = alpha;
        double Energy_candidate;
        SCFEnergy SCF_Energy_candidate;
        OrbitalVector dPhi_n;
        while (true) {
            F.clear();
            // Retraction to Stiefel is Lowdin based:
            Phi_n = orbital::add(1.0, Phi_backup, alpha_trial, direction);
            // Orthonormalization updates F_mat as a side effect?!
            orbital::orthonormalize(orb_prec, Phi_n, F_mat);
            // Compute Fock matrix and energy
            if (F.getReactionOperator() != nullptr) F.getReactionOperator()->updateMOResidual(err_t);
            F.setup(orb_prec);
            F_mat = F(Phi_n, Phi_n);
            SCF_Energy_candidate = F.trace(Phi_n, nucs);
            Energy_candidate = SCF_Energy_candidate.getTotalEnergy();
            println(0, "Candidate energy: " << Energy_candidate);

            dPhi_n = orbital::add(1.0, Phi_n, -1.0, Phi_backup);
            errors = orbital::get_norms(dPhi_n);
            dPhi_n.clear();
            err_o = errors.maxCoeff();
            if (checkConvergence(err_o, 0.0))
            {
                println(0, "Line search step is negligible; accepting.");
                break;
            }

            // Directional Armijo condition:
            if (Energy_candidate <= Energy + armijo_parameter * alpha_trial * descent_directional_product) {
                // Accept step
                println(0, "update: " << err_o << " (step size = " << alpha_trial << ")");
                break;
            }
            alpha_trial *= beta;
            if (alpha_trial < alpha_min) {
                println(0, "Warning: step size too small, stopping search.");
                break;
            }
            rejectness_count += 1;
        }

        // Step-size growth safeguard
        if (Energy_candidate <= Energy + 0.7 * alpha_trial * descent_directional_product) {
            alpha = std::min(alpha_trial * gamma, alpha_max);   // grow step size
        } else {
            alpha = alpha_trial;                          // keep shrunk step size
        }
    

        // Compute total update norm and collect convergence data
        err_t = errors.norm();
        json_cycle["mo_residual"] = err_t;
        this->error.push_back(err_t);

        // Energy_candidate < Energy, unless convergence is reached
        if (Energy_candidate >= Energy) {
            println(0, "Energy cannot be decreased more in the backtracking search.");
            converged = true;
            Phi_n = Phi_backup;
        }
        else {
            E_n = SCF_Energy_candidate;
        }
        this->energy.push_back(E_n);
        this->property.push_back(E_n.getTotalEnergy());
        Phi_backup.clear();

        auto err_p = calcPropertyError();
        if (not converged) converged = checkConvergence(err_o, err_p);

        json_cycle["energy_terms"] = E_n.json();
        json_cycle["energy_total"] = E_n.getTotalEnergy();
        json_cycle["energy_update"] = err_p;

        mrcpp::print::separator(0, '-');

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

    // Print energies, gradients, properties for debugging
    println(0, "The amount of Armijo rejections: " << rejectness_count);
    println(0, "norm(grad_E):");
    for (std::size_t i = 0; i < grad_E_array.size(); ++i)
        println(0, "norm(grad_E)[" << i << "] = " << grad_E_array[i]);
    

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

/** 
 * @brief Determine new occupation vector according to MOM/IMOM procedure
 * @param Phi_n: orbitals of current iteration n.
 * @param Phi_mom: orbitals of last iteration n-1 (MOM) or first iteration (IMOM).
 * 
 * According to MOM/IMOM procedure the occupation numbers for the current iteration get
 * determined based on the overlap with the orbitals of an earlier iteration of the SCF procedure.
 */
DoubleVector GroundStateSolver::getNewOccupations(OrbitalVector &Phi_n, OrbitalVector &Phi_mom) {
    DoubleMatrix overlap = orbital::calc_overlap_matrix(Phi_mom, Phi_n).real();
    DoubleVector occup = orbital::get_occupations(Phi_mom); // get occupation numbers of the orbitals of the first iteration
    double occ1 = occup(0);
    DoubleVector occNew = DoubleVector::Constant(occup.size(), occ1);

    // create vector which contains the positions of the second occupation number
    DoubleVector currOcc = DoubleVector::Zero(occup.size());
    unsigned int nCurrOcc = 0;
    double occ2 = 0.0;
    for (unsigned int i = 1; i < occup.size(); i++) {
        if (occup(i) != occ1) {
            occ2 = occup(i);
            currOcc(i) = 1.0;
            nCurrOcc++;
        }
    }

    // only consider overlap with orbitals with the second occupation number
    DoubleMatrix occOverlap = currOcc.asDiagonal() * overlap;
    DoubleVector p = occOverlap.colwise().norm();

    // debug print section
    print_utils::matrix(3, "MOM overlap matrix", overlap, 2);
    print_utils::vector(3, "MOM total overlap", p, 2);

    // sort by highest overlap
    std::vector<std::pair<double, unsigned int>> sortme;
    for (unsigned int q = 0; q < p.size(); q++)
        sortme.push_back(std::pair<double, unsigned int>(p(q), q));
    std::stable_sort(sortme.begin(), sortme.end());
    std::reverse(sortme.begin(), sortme.end());

    // assign the second occupation number to orbitals with highest overlap
    for (unsigned int q = 0; q < nCurrOcc; q++)
        occNew(sortme[q].second) = occ2;
    return occNew;
}

} // namespace mrchem
