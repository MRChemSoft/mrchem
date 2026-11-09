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

////////////////
// Debugging
///////////////
#include "qmoperators/one_electron/H_B_dip.h"
#include "qmoperators/one_electron/H_BB_dia.h"

namespace mrchem {

static void testMagneticOperators(double prec,
                                   std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                                   const mrcpp::Coord<3> &o) {
    mrcpp::print::separator(0, '=');
    println(0, "  Testing H_B_dip[2] and H_BB_dia[2][2] on Gaussians");
    mrcpp::print::separator(0, '-');

    const double alpha = 1.0;

    // Helper: norm of a FunctionTree (getSquareNorm() returns ||f||^2)
    auto treeNorm = [](const mrcpp::FunctionTree<3> &t) {
        return std::sqrt(t.getSquareNorm());
    };

    // =========================================================
    // Test 1: H_BB_dia[2][2] on phi = exp(-alpha*r^2)
    // Expected (real): (1/4)*((x-ox)^2+(y-oy)^2)*exp(-alpha*r^2)
    // =========================================================
    {
        Orbital phi(SPIN::Alpha);
        phi.alloc(NUMBER::Real);
        std::function<double(const mrcpp::Coord<3> &)> f_phi =
            [alpha](const mrcpp::Coord<3> &r) -> double {
            return std::exp(-alpha * (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
        };
        mrcpp::project<3>(prec, phi.real(), f_phi);  // <-- explicit template arg

        H_BB_dia h_dia(o);
        Orbital result = h_dia[2][2](phi);

        auto &MRA = phi.real().getMRA();
        mrcpp::FunctionTree<3> expected(MRA);
        std::function<double(const mrcpp::Coord<3> &)> f_exp =
            [alpha, o](const mrcpp::Coord<3> &r) -> double {
            double dx = r[0]-o[0], dy = r[1]-o[1];
            return 0.25*(dx*dx+dy*dy)*std::exp(-alpha*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
        };
        mrcpp::project<3>(prec, expected, f_exp);

        double nrm = result.hasReal() ? treeNorm(result.real()) : 0.0;
        double err = -1.0;
        if (result.hasReal()) {
            mrcpp::FunctionTree<3> diff(MRA);
            mrcpp::add(prec, diff, 1.0, result.real(), -1.0, expected);
            err = treeNorm(diff);
        }
        println(0, "  H_BB_dia[2][2] on exp(-r^2):");
        println(0, "    ||result_real||  = " << nrm << "  (expected ~0.741)");
        println(0, "    ||imag||         = " << (result.hasImag() ? treeNorm(result.imag()) : 0.0) << "  (expected: 0)");
        println(0, "    ||error||        = " << err);
        println(0, "    relative err     = " << (nrm > 0 ? err/nrm : -1.0));

        phi.free(NUMBER::Total);
        result.free(NUMBER::Total);
    }

    // =========================================================
    // Test 2: H_B_dip[2] on phi = x*exp(-alpha*r^2)
    // Expected: real = 0,  imag = +(1/2)*y*exp(-alpha*r^2)
    // =========================================================
    {
        Orbital phi(SPIN::Alpha);
        phi.alloc(NUMBER::Real);
        std::function<double(const mrcpp::Coord<3> &)> f_phi =
            [alpha](const mrcpp::Coord<3> &r) -> double {
            return r[0]*std::exp(-alpha*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
        };
        mrcpp::project<3>(prec, phi.real(), f_phi);

        H_B_dip h_dip(D, o);
        Orbital result = h_dip[2](phi);

        auto &MRA = phi.real().getMRA();
        mrcpp::FunctionTree<3> expected_imag(MRA);
        std::function<double(const mrcpp::Coord<3> &)> f_exp =
            [alpha](const mrcpp::Coord<3> &r) -> double {
            return 0.5*r[1]*std::exp(-alpha*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
        };
        mrcpp::project<3>(prec, expected_imag, f_exp);

        double real_nrm = result.hasReal() ? treeNorm(result.real()) : 0.0;
        double imag_nrm = result.hasImag() ? treeNorm(result.imag()) : 0.0;
        double err = -1.0;
        if (result.hasImag()) {
            mrcpp::FunctionTree<3> diff(MRA);
            mrcpp::add(prec, diff, 1.0, result.imag(), -1.0, expected_imag);
            err = treeNorm(diff);
        }
        println(0, "  H_B_dip[2] on x*exp(-r^2):");
        println(0, "    ||real part||    = " << real_nrm << "  (expected: 0)");
        println(0, "    ||imag part||    = " << imag_nrm << "  (expected ~0.370)");
        println(0, "    ||imag error||   = " << err);

        phi.free(NUMBER::Total);
        result.free(NUMBER::Total);
    }

    // =========================================================
    // Test 3: H_B_dip[2] on phi = exp(-alpha*r^2), gauge o=(0,1,0)
    // Expected: real = 0,  imag = alpha*x*exp(-alpha*r^2)
    // =========================================================
    {
        mrcpp::Coord<3> o_offset = {0.0, 1.0, 0.0};
        Orbital phi(SPIN::Alpha);
        phi.alloc(NUMBER::Real);
        std::function<double(const mrcpp::Coord<3> &)> f_phi =
            [alpha](const mrcpp::Coord<3> &r) -> double {
            return std::exp(-alpha*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
        };
        mrcpp::project<3>(prec, phi.real(), f_phi);

        H_B_dip h_dip(D, o_offset);
        Orbital result = h_dip[2](phi);

        auto &MRA = phi.real().getMRA();
        mrcpp::FunctionTree<3> expected_imag(MRA);
        std::function<double(const mrcpp::Coord<3> &)> f_exp =
            [alpha](const mrcpp::Coord<3> &r) -> double {
            return alpha*r[0]*std::exp(-alpha*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
        };
        mrcpp::project<3>(prec, expected_imag, f_exp);

        double real_nrm = result.hasReal() ? treeNorm(result.real()) : 0.0;
        double imag_nrm = result.hasImag() ? treeNorm(result.imag()) : 0.0;
        double err = -1.0;
        if (result.hasImag()) {
            mrcpp::FunctionTree<3> diff(MRA);
            mrcpp::add(prec, diff, 1.0, result.imag(), -1.0, expected_imag);
            err = treeNorm(diff);
        }
        println(0, "  H_B_dip[2] on exp(-r^2), gauge=(0,1,0):");
        println(0, "    ||real part||    = " << real_nrm << "  (expected: 0)");
        println(0, "    ||imag part||    = " << imag_nrm << "  (expected ~0.741)");
        println(0, "    ||imag error||   = " << err);

        phi.free(NUMBER::Total);
        result.free(NUMBER::Total);
    }

    mrcpp::print::separator(0, '=');
}

} // namespace mrchem
////////////////
////////////////

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
    print_utils::text(0, "XC library         ", (this->libxc) ? "LibXC" : "XCFun");
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

    // Convert to complex if magnetic perturbation is active
    if (not F.perturbation().name().empty()) {
        Phi_n = orbital::CopyToComplex(Phi_n);
    }

    //////////////////////////////
    // Debugging: test magnetic operators
    if (not F.perturbation().name().empty()) {
        auto D_test = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
        mrcpp::Coord<3> o_test = {0.0, 0.0, 0.0};
        testMagneticOperators(this->orbPrec[0], D_test, o_test);
    }
    //////////////////////////////

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

        // Diagonalize F_mat
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigensolver(F_mat);
        if (eigensolver.info() != Eigen::Success) {
            MSG_ABORT("Eigen-decomposition of F_mat failed");
        }

        Eigen::VectorXd sigma_F = eigensolver.eigenvalues();
        ComplexMatrix U_F = eigensolver.eigenvectors();

        // Init Helmholtz operator
        HelmholtzVector H(helm_prec, sigma_F);
        //ComplexMatrix L_mat = H.getLambdaMatrix();

        // Apply Helmholtz operator
        //OrbitalVector Psi = F.buildHelmholtzArgument(orb_prec, Phi_n, F_mat, L_mat);
        OrbitalVector Psi = F.potential()(Phi_n);
        Psi = orbital::rotate(Psi, U_F.adjoint(), orb_prec);
        OrbitalVector Phi_np1 = H(Psi);
        Psi.clear();
        Phi_np1 = orbital::rotate(Phi_np1, U_F, orb_prec);


        // Born-Neumann second term: paramagnetic  R_0 * W_para * R_0 * source
        // Requires complex orbitals — enabled after CopyToComplex above
        if (not F.perturbation().name().empty()) {
            OrbitalVector W_para_Phi = F.perturbation()(Phi_np1);

            // H_B_dip stores result in real part (isImag not honored by CompFunction).
            // Multiply by i: orbital::add with ComplexDouble coefficient correctly
            // routes  real→imag and imag→−real, handling both components.
            OrbitalVector zeros = orbital::param_copy(W_para_Phi);  // right structure, no trees
            W_para_Phi = orbital::add(ComplexDouble(0.0, 1.0), W_para_Phi,
                                    ComplexDouble(0.0, 0.0), zeros, orb_prec);
            zeros.clear();

            W_para_Phi = orbital::rotate(W_para_Phi, U_F.adjoint(), orb_prec);
            OrbitalVector corr = H(W_para_Phi);
            W_para_Phi.clear();
            corr = orbital::rotate(corr, U_F, orb_prec);
            for (int i = 0; i < static_cast<int>(Phi_np1.size()); i++)
                if (mrcpp::mpi::my_func(Phi_np1[i])) Phi_np1[i].add(1.0, corr[i]);
            corr.clear();
        }

        // Born-Neumann second term: Phi += R_0 * W_dia * R_0 * source  (real, diamagnetic only)
        // The paramagnetic term (W_para) maps real->imaginary and requires complex orbital
        // support in MRCPP; it is skipped until that is available.
        if (not F.diaField().name().empty()) {
            OrbitalVector W_dia_Phi = F.diaField()(Phi_np1);
            W_dia_Phi = orbital::rotate(W_dia_Phi, U_F.adjoint(), orb_prec);
            OrbitalVector corr      = H(W_dia_Phi);
            W_dia_Phi.clear();
            corr = orbital::rotate(corr, U_F, orb_prec);
            for (int i = 0; i < static_cast<int>(Phi_np1.size()); i++) {
                if (mrcpp::mpi::my_func(Phi_np1[i]))
                    Phi_np1[i].add(1.0, corr[i]);
            }
            corr.clear();
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

        // MOM / IMOM: get the new occupation vector for the current scf iteration
        OrbitalVector &Phi_mom = mol.getOrbitalsMom();
        if (Phi_mom.size() > 0) {
            bool restricted = (orbital::size_doubly(Phi_n) != 0);
            if (restricted) {
                DoubleVector occNew = getNewOccupations(Phi_n, Phi_mom);
                orbital::set_occupations(Phi_n, occNew);
                mol.calculateOrbitalPositions();
                if (plevel >= 2)
                    mol.printOrbitalPositions();
            }
            else {
                // in case of unrestricted calculation, get the new occupation for alpha and beta spins independently
                OrbitalVector Phi_n_beta = orbital::deep_copy(Phi_n);
                OrbitalVector Phi_mom_beta = orbital::deep_copy(Phi_mom);
                OrbitalVector Phi_n_a = orbital::disjoin(Phi_n_beta, SPIN::Alpha);
                OrbitalVector Phi_mom_a = orbital::disjoin(Phi_mom_beta, SPIN::Alpha);
                DoubleVector occAlpha = getNewOccupations(Phi_n_a, Phi_mom_a);
                DoubleVector occBeta = getNewOccupations(Phi_n_beta, Phi_mom_beta);
                DoubleVector occNew(occAlpha.size() + occBeta.size());
                occNew << occAlpha, occBeta;
                orbital::set_occupations(Phi_n, occNew);
                mol.calculateOrbitalPositions();
                if (plevel >= 2)
                    mol.printOrbitalPositions();
            }
        }
        // MOM: save orbitals of current iteration for next iteration of the SCF procedure
        if (deltaSCFMethod == "MOM")
            Phi_mom = orbital::deep_copy(Phi_n);

        orbital::orthonormalize(orb_prec, Phi_n, F_mat);

        // Compute Fock matrix and energy
        if (F.getReactionOperator() != nullptr) F.getReactionOperator()->updateMOResidual(err_t);
        F.setup(orb_prec);
        F_mat = F(Phi_n, Phi_n);
        E_n = F.trace(Phi_n, nucs);

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
