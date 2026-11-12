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

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "LagrangianSolver.h"

#include "HelmholtzVector.h"
#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/FockBuilder.h"
#include "utils/print_utils.h"
#include <iostream>

using mrcpp::Printer;
using mrcpp::Timer;
using nlohmann::json;
using PoissonOperator = mrcpp::PoissonOperator;

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

namespace mrchem {

namespace {
/** @brief Convert all orbitals to complex if at least one of them is complex (or if forced)
 *
 * Complex coefficients in the orbital update can turn some orbitals complex and
 * leave others real, while MRCPP routines (e.g. rotate, in-place add) expect a uniform type.
 */
void make_uniform_type(OrbitalVector &Phi, bool force_complex = false) {
    bool any_complex = force_complex;
    for (auto &phi : Phi)
        if (phi.iscomplex()) any_complex = true;
    if (!any_complex) return;

    for (auto &phi : Phi) {
        if (phi.iscomplex()) continue;
        Orbital phi_c;
        mrcpp::CopyToComplex(phi_c, phi);
        phi = phi_c;
    }
}
} // namespace

// TODO: read all the settings in input
LagrangianSolver::LagrangianSolver(double prec) {
    this->prec = prec;
    this->threshold = 1e-10;
    this->P_p = std::make_shared<PoissonOperator>(*MRA, this->prec);
    this->history = 0;
}

void LagrangianSolver::set_orbitals(OrbitalVector Phi_n) {
    this->orbitals = std::make_shared<OrbitalVector>(Phi_n);
}

/** @brief Switch to complex orbitals (and complex integrals and RDMs in the external solver) */
void LagrangianSolver::set_complex_orbitals(ChemTensorSolver &S) {
    this->real_orbitals = false;
    make_uniform_type(*this->orbitals, true);
    S.set_real(false);
}

/** @brief Run Lagrangian orbital optimization
 *
 * @param mol: Molecule to optimize
 * @param F: Fock Builder
 * @param S: External solver to compute 1-2RDMs
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm performs
 * a Lagrangian optimization of the orbitals, where the constraint is the orthonormal
 * condition of the orbitals. The driver is used to compute the 1-2RDMs.
 */
json LagrangianSolver::optimize(Molecule &mol, FockBuilder &F, ChemTensorSolver &S) {
    // printParameters("Optimize ground state orbitals");

    Timer t_tot;
    json json_out;

    // Nuclei nucs = mol.getNuclei();

    // initial orbitals
    set_orbitals(mol.getOrbitals());
    // real orbitals (cheaper) whenever possible: complex if requested or if any initial orbital is complex.
    // The calculation switches to complex later if the integrals turn out complex
    bool any_complex = false;
    for (auto &phi : *this->orbitals)
        if (phi.iscomplex()) any_complex = true;
    this->real_orbitals = true;
    S.set_real(true);
    if (this->force_complex or any_complex) set_complex_orbitals(S);
    print_utils::text(0, "Orbitals", this->real_orbitals ? "real" : "complex");
    DoubleVector errors = DoubleVector::Ones(orbitals->size());

    if (this->localize) {
        ComplexMatrix tmpfock(orbitals->size(), orbitals->size());
        orbital::localize(this->prec, *this->orbitals, tmpfock);
    }

    std::shared_ptr<KAIN> kain;

    if (this->history) {
        auto scaling = std::vector<double>(orbitals->size(), 1.0);
        kain = std::make_shared<KAIN>(this->history, 0, false, scaling);
    }

    bool converged = false;
    json_out["cycles"] = {};
    for (int nIter = 0; nIter < this->maxIter; nIter++) {
        json json_cycle;
        std::stringstream o_header;
        o_header << "Lagrangian cycle " << nIter;
        mrcpp::print::header(1, o_header.str(), 0, '#');
        Timer t_scf;

        // set integrals, run DMRG and calculate RDMS
        S.set_integrals(*this->orbitals);
        // the external solver switches to the complex mode if the integrals are complex
        if (this->real_orbitals and not S.is_real()) set_complex_orbitals(S);
        S.optimize();

        // check if converged (before the update, so that orbitals and last energy match)
        // NOTE: no orbital threshold, convergence only on the energy
        this->property.push_back(S.get_energy());
        double err_p = getUpdate(this->property, this->property.size(), true);
        converged = (nIter > 0) and checkConvergence(0.0, err_p);

        json_cycle["energy_total"] = S.get_energy();
        json_cycle["energy_update"] = err_p;
        auto pprec = Printer::getPrecision();
        print_utils::scalar(0, "Total energy", S.get_energy(), "(au)", 2 * pprec);
        print_utils::scalar(0, "Energy difference", err_p, "(au)", pprec, true);

        // update orbitals (not needed if converged or in the last iteration)
        if (not converged and nIter < this->maxIter - 1) {
            errors = orbital_update(F, S, kain);
            make_uniform_type(*this->orbitals);

            // reinforce orthonormalization immediately (Loewdin: Phi = S^(-1/2) Phi)
            ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(*this->orbitals);
            if (this->real_orbitals) S_m12 = S_m12.real().cast<ComplexDouble>();
            mrcpp::rotate(*this->orbitals, S_m12, this->prec);

            if (this->real_orbitals)
                for (auto &phi : *this->orbitals)
                    if (phi.iscomplex()) MSG_ABORT("Complex orbital in the real mode");
        }

        printOrbitals(S.get_lagrange_multipliers()->real().diagonal(), S.get_one_rdm()->real().diagonal(), errors, *this->orbitals);

        t_scf.stop();
        json_cycle["wall_time"] = t_scf.elapsed();
        mrcpp::print::footer(1, t_scf, 2, '#');

        json_out["cycles"].push_back(json_cycle);
        if (converged) break;
    }

    F.clear();

    auto w0 = Printer::getWidth() - 1;
    mrcpp::print::separator(0, '-');
    if (converged) {
        println(0, std::string((w0 - 45) / 2, ' ') << "Lagrangian optimization converged in " << this->property.size() - 1 << " iterations!");
    } else {
        println(0, std::string((w0 - 40) / 2, ' ') << "Lagrangian optimization did NOT converge!!!");
    }
    mrcpp::print::separator(0, '=', 2);
    reset();

    json_out["wall_time"] = t_tot.elapsed();
    json_out["converged"] = converged;
    return json_out;
}

DoubleVector LagrangianSolver::orbital_update(FockBuilder &F, ChemTensorSolver &S, std::shared_ptr<KAIN> kain) {
    mrcpp::print::header(3, "Computing orbital update", 0, '-');
    Timer t_tot, t_1e(false), t_2e(false);

    const int L = this->orbitals->size();

    // calculate Lagrange multipliersin the old basis
    // (otherwise you need also two-boyd integral and rdm in the new basis)
    S.calculate_lagrange_multipliers();

    // diagonalize 1rdm. change basis of one-body integral, multipliers and orbitals
    // (two-body integral and rdms not necessary)
    S.diagonalize_1rdm();
    std::shared_ptr<ComplexMatrix> basis_change = S.get_basis_change();

    // initialize Helmholtz operator
    S.calculate_helmholtz_coefficients();
    std::shared_ptr<DoubleVector> ptr_helm_coeff = S.get_helmholtz_coefficients();
    // as in python: positive coefficients (mu imaginary) are not allowed. The check is done here,
    // since HelmholtzVector silently replaces them with -0.5 (needed by the standard SCF solver)
    for (int m = 0; m < L; m++) {
        double x_m = (*ptr_helm_coeff)(m);
        if (x_m > this->threshold) {
            std::stringstream msg;
            msg << "Positive Helmholtz coefficient for orbital " << m << ": " << x_m;
            MSG_ABORT(msg.str());
        }
    }
    HelmholtzVector H(this->prec, *ptr_helm_coeff);

    // new orbitals
    OrbitalVector new_Phi;

    t_1e.resume();
    orbital_update_one_body(F, S, new_Phi);
    t_1e.stop();

    t_2e.resume();
    orbital_update_two_body(F, S, new_Phi);
    t_2e.stop();

    t_tot.stop();

    mrcpp::print::time(3, "Time 1e terms", t_1e);
    mrcpp::print::time(3, "Time 2e terms", t_2e);
    mrcpp::print::time(3, "Time updating orbitals", t_tot);
    mrcpp::print::separator(3, '-');

    // apply Helmholtz kernel and update orbitals
    OrbitalVector new_orbitals = H(new_Phi);
    mrcpp::rotate(new_orbitals, *basis_change); // rotate back to old basis
    OrbitalVector diff_orbitals = orbital::add(1.0, new_orbitals, -1.0, *this->orbitals);

    if (kain) { kain->accelerate(1e-10, *this->orbitals, diff_orbitals); }

    new_orbitals = orbital::add(1.0, *this->orbitals, 1.0, diff_orbitals);

    DoubleVector errors = orbital::get_norms(diff_orbitals);
    (*this->orbitals) = new_orbitals;
    return errors;
}

void LagrangianSolver::orbital_update_one_body(FockBuilder &F, ChemTensorSolver &S, OrbitalVector &new_Phi) {
    const int L = this->orbitals->size();
    NuclearOperator &V = *(F.getNuclearOperator());
    ComplexMatrix &lag_coeff = *(S.get_lagrange_multipliers());
    ComplexMatrix &one_rdm = *(S.get_one_rdm());
    ComplexMatrix &basis_change = *(S.get_basis_change());

    OrbitalVector V_phi = V(*this->orbitals);
    mrcpp::rotate(V_phi, basis_change.adjoint());

    ComplexMatrix rot_coef(L, L);
    for (size_t i = 0; i < L; i++) {
        ComplexDouble denom = one_rdm(i, i);
        for (size_t j = 0; j < L; j++) { rot_coef(i, j) = coeff(i != j ? -lag_coeff(i, j) / denom : 0.0); }
    }

    ComplexMatrix net_rot_matrix = rot_coef * basis_change;

    OrbitalVector rotated = orbital::rotate(*this->orbitals, net_rot_matrix.adjoint());
    new_Phi = orbital::add(1.0, V_phi, 1.0, rotated);
}

void LagrangianSolver::orbital_update_two_body(FockBuilder &F, ChemTensorSolver &S, OrbitalVector &new_Phi) {
    const int L = this->orbitals->size();
    auto &one_rdm = *(S.get_one_rdm());
    auto &two_rdm = *(S.get_two_rdm());
    auto &lag_coeff = *(S.get_lagrange_multipliers());
    auto &basis_change = *(S.get_basis_change());
    auto &P = *this->P_p;
    auto &orbitals = *this->orbitals;

    // precompute U_dE2(m,j,k,l) = sum_i U(m,i) * two_rdm(i,j,k,l)
    ComplexTensorR4 U_dE2(L, L, L, L);
    U_dE2.setZero();
    for (int m = 0; m < L; m++)
        for (int i = 0; i < L; i++)
            for (int j = 0; j < L; j++)
                for (int k = 0; k < L; k++)
                    for (int l = 0; l < L; l++) U_dE2(m, j, k, l) += basis_change(m, i) * two_rdm(i, j, k, l);

    // two-body term: m on the new basis (through U_dE2), j, k, l on the old basis (like two_rdm)
    for (int j = 0; j < L; j++) {
        std::cout << "Orbital " << j + 1 << " out of " << L << std::endl;
        for (int l = j; l < L; l++) { // Triangular loop l >= j
            // g_jl = poisson(Phi[j]^* * Phi[l]) (the MRCPP kernel is already 1/|r-r'|, no 4*pi as in python)
            Orbital phi_jl, g_jl;
            mrcpp::multiply(phi_jl, orbitals[j], orbitals[l], this->prec, false, false, true);

            if (phi_jl.norm() < this->prec * 1.0e-5) {
                // std::cout << "Skipping orbital pair: " << j << ", " << l << std::endl;
                continue;
            }

            mrcpp::apply(this->prec, g_jl, P, phi_jl);

            for (int k = 0; k < L; k++) {
                // state = g_jl * Phi[k]
                Orbital state;
                mrcpp::multiply(state, g_jl, orbitals[k], this->prec);

                for (int m = 0; m < L; m++) {
                    ComplexDouble denom;
                    if (std::abs(one_rdm(m, m)) > threshold)
                        denom = one_rdm(m, m);
                    else if (std::abs(lag_coeff(m, m)) > threshold)
                        denom = lag_coeff(m, m);
                    else
                        continue; // both below threshold, skip

                    ComplexDouble u = U_dE2(m, j, k, l);
                    if (l != j) u += U_dE2(m, l, k, j);

                    new_Phi[m].add(coeff(2.0 * u / denom), state);
                }
            }
        }
    }
}

/** @brief Pretty printing of orbitals with energies and occupancies
 *
 * @param epsilon: orbital energies
 * @param Phi: orbital vector
 */
void LagrangianSolver::printOrbitals(const DoubleVector &energies, const DoubleVector &occupancies, const DoubleVector &errors, OrbitalVector &Phi, bool print_head) const {
    int pprec = Printer::getPrecision();
    int w0 = (Printer::getWidth() - 1);
    int w1 = 5;
    int w2 = 7;
    int w3 = 19;
    int w4 = 19;
    int w5 = 17;
    int w6 = w0 - w1 - w2 - w3 - w4 - w5;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w2) << "Nodes";
    o_head << std::setw(w3) << "Energy";
    o_head << std::setw(w4) << "Occupancy";
    o_head << std::setw(w5) << "Residual";
    o_head << std::setw(w6) << "Done";

    if (print_head) mrcpp::print::separator(2, '=');
    if (print_head) println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    for (size_t i = 0; i < Phi.size(); i++) {
        bool conv_i = (errors(i) < this->orbThrs) or (this->orbThrs < 0.0);
        std::stringstream o_row;
        o_row << std::setw(w1) << i;
        o_row << std::setw(w2) << Phi[i].getNNodes();
        o_row << std::setw(w3) << std::setprecision(2 * pprec) << std::fixed << energies(i) * 0.5;
        o_row << std::setw(w4) << std::setprecision(2 * pprec) << std::fixed << occupancies(i);
        o_row << std::setw(w5) << std::setprecision(pprec) << std::scientific << errors(i);
        o_row << std::setw(w6) << conv_i;
        println(2, o_row.str());
        if (Phi[i].hasReal()) println(5, Phi[i].real());
    }
}

} // namespace mrchem