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

        // Apply Helmholtz operator
        OrbitalVector Psi = F.buildHelmholtzArgument(orb_prec, Phi_n, F_mat, L_mat);
        OrbitalVector Phi_np1 = H(Psi);
        Psi.clear();
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
