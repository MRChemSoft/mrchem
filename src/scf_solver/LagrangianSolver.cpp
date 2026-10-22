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
#include <MRCPP/MWOperators>

#include "LagrangianSolver.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/two_electron/FockBuilder.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "HelmholtzVector.h"
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
LagrangianSolver::LagrangianSolver(double prec){
    this->prec = prec;
    this->threshold = 1e-10;
    this->P_p = std::make_shared<PoissonOperator>(*MRA, this->prec);
}

void LagrangianSolver::set_orbitals(OrbitalVector Phi_n){
    this->orbitals = std::make_shared<OrbitalVector>(Phi_n);
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
    //printParameters("Optimize ground state orbitals");

    Timer t_tot;
    json json_out;

    //Nuclei nucs = mol.getNuclei();

    // initial orbitals
    set_orbitals(mol.getOrbitals());
    // work with complex orbitals throughout: the RDMs from the external solver are complex,
    // and the in-place add in the orbital update needs orbitals and coefficients of the same type
    // TODO: distinguish between real and complex cases (real orbitals are cheaper)
    make_uniform_type(*this->orbitals, true);

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
            orbital_update(F, S);
            make_uniform_type(*this->orbitals);

            // reinforce orthonormalization immediately (Loewdin: Phi = S^(-1/2) Phi)
            ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(*this->orbitals);
            mrcpp::rotate(*this->orbitals, S_m12, this->prec);
        }

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

void LagrangianSolver::orbital_update(FockBuilder &F, ChemTensorSolver &S){
    const int L = this->orbitals->size();

    // calculate Lagrange multipliersin the old basis
    // (otherwise you need also two-boyd integral and rdm in the new basis)
    S.calculate_lagrange_multipliers();

    // diagonalize 1rdm. change basis of one-body integral, multipliers and orbitals 
    // (two-body integral and rdms not necessary)
    S.diagonalize_1rdm();
    std::shared_ptr<ComplexMatrix> basis_change = S.get_basis_change();
    // keep the old orbitals: the two-body term uses them together with the old-basis two-body rdm
    OrbitalVector old_Phi = *this->orbitals;
    orbital_basis_change(basis_change);

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
    OrbitalVector new_Phi(L);

    orbital_update_one_body(F, S, new_Phi);
    orbital_update_two_body(F, S, old_Phi, new_Phi);
    
    // apply Helmholtz kernel and update orbitals
    (*this->orbitals) = H(new_Phi);
}

void LagrangianSolver::orbital_basis_change(std::shared_ptr<ComplexMatrix> basis_change){
    const int L = this->orbitals->size();
    OrbitalVector new_Phi;

    for (int i = 0; i < L; i++) {
        // phi_i = sum_j U(i,j) * Phi[j]
        std::vector<ComplexDouble> coeffs(L);
        for (int j = 0; j < L; j++)
            coeffs[j] = (*basis_change)(i, j);
        Orbital phi_i;
        mrcpp::linear_combination(phi_i, coeffs, *this->orbitals, this->prec);
        new_Phi.push_back(phi_i);
    }

    *this->orbitals = new_Phi;
}

void LagrangianSolver::orbital_update_one_body(FockBuilder &F, ChemTensorSolver &S, OrbitalVector &new_Phi){
    const int L = this->orbitals->size();
    NuclearOperator &V = *(F.getNuclearOperator());
    ComplexMatrix &lag_coeff = *(S.get_lagrange_multipliers());
    ComplexMatrix &one_rdm = *(S.get_one_rdm());

    for(auto m=0; m<L; m++){
        if (std::abs(one_rdm(m, m)) > this->threshold) {
            // nuclear potential term (on new basis)
            new_Phi[m] = V((*this->orbitals)[m]);
            //Orbital tmp = V((*this->orbitals)[m]);
            //mrcpp::add(new_Phi[m], 1.0, new_Phi[m], 1.0, tmp, this->prec);
            for (auto j = 0; j < L; j++) {
                // multipliers' term (on new basis)
                if (j != m)
                    new_Phi[m].add(-lag_coeff(m, j) / one_rdm(m, m), (*this->orbitals)[j]);
            }
        } else if (std::abs(lag_coeff(m, m)) > this->threshold) {
            std::cerr << "WARNING: lambda_m below threshold, but not helm_coeff'_mm !" << std::endl;
            // start from zero
            new_Phi[m] = (*this->orbitals)[m].paramCopy(true);
            for (auto j = 0; j < L; j++) {
                if (j != m)
                    new_Phi[m].add(-lag_coeff(m, j) / lag_coeff(m, m), (*this->orbitals)[j]);
            }
        } else {
            std::cerr << "WARNING: Both lambda_m and helm_coeff'_mm below threshold" << std::endl;
            new_Phi[m] = (*this->orbitals)[m];
        }
    }
}

void LagrangianSolver::orbital_update_two_body(FockBuilder &F, ChemTensorSolver &S, OrbitalVector &old_Phi, OrbitalVector &new_Phi){
    const int L = this->orbitals->size();
    auto &one_rdm = *(S.get_one_rdm());
    auto &two_rdm = *(S.get_two_rdm());
    auto &lag_coeff = *(S.get_lagrange_multipliers());
    auto &basis_change = *(S.get_basis_change());
    auto &P = *this->P_p;

    // precompute U_dE2(m,j,k,l) = sum_i U(m,i) * two_rdm(i,j,k,l)
    ComplexTensorR4 U_dE2(L, L, L, L);
    U_dE2.setZero();
    for (int m = 0; m < L; m++)
        for (int i = 0; i < L; i++)
            for (int j = 0; j < L; j++)
                for (int k = 0; k < L; k++)
                    for (int l = 0; l < L; l++)
                        U_dE2(m, j, k, l) += basis_change(m, i) * two_rdm(i, j, k, l);

    // two-body term: m on the new basis (through U_dE2), j, k, l on the old basis (like two_rdm)
    for (int j = 0; j < L; j++) {
        std::cout << "Orbital " << j+1 << " out of " << L << std::endl;
        for (int l = 0; l < L; l++) {
            // g_jl = poisson(Phi[j]^* * Phi[l]) (the MRCPP kernel is already 1/|r-r'|, no 4*pi as in python)
            Orbital phi_jl, g_jl;
            mrcpp::multiply(phi_jl, old_Phi[j], old_Phi[l], this->prec, false, false, true);
            mrcpp::apply(this->prec, g_jl, P, phi_jl);

            for (int k = 0; k < L; k++) {
                // state = g_jl * Phi[k]
                Orbital state;
                mrcpp::multiply(state, g_jl, old_Phi[k], this->prec);

                for (int m = 0; m < L; m++) {
                    ComplexDouble denom;
                    if (std::abs(one_rdm(m, m)) > threshold)
                        denom = one_rdm(m, m);
                    else if (std::abs(lag_coeff(m, m)) > threshold)
                        denom = lag_coeff(m, m);
                    else
                        continue; // both below threshold, skip

                    ComplexDouble coeff = 2.0 * U_dE2(m, j, k, l) / denom;
                    new_Phi[m].add(coeff, state);
                }
            }
        }
    }
}




} // namespace mrchem