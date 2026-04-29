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
#include <iostream>

using mrcpp::Printer;
//using mrcpp::Timer;
using nlohmann::json;
using PoissonOperator = mrcpp::PoissonOperator;

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;


namespace mrchem {

// TODO: read all the settings in input
LagrangianSolver::LagrangianSolver(){
    this->P_p = std::make_shared<PoissonOperator>(*MRA, this->prec);
    this->nIter = 10;
    this->scf_tol = 1e-3;
    this->prec = 1e-3;
    this->threshold = 1e-10;
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

    json json_out;

    //Nuclei nucs = mol.getNuclei();

    // initial orbitals
    set_orbitals(mol.getOrbitals());
    double prec = 1e-3;

    bool converged = false;
    
    for(int i=0; i<this->nIter && !converged; i++){
        // set integrals, run DMRG and calculate RDMS
        S.set_integrals(*this->orbitals);
        S.optimize();
        this->energy.push_back(S.get_energy());
        std::cout << "Energy:" << S.get_energy() << std::endl;

        // update orbitals
        orbital_update(F, S);

        // check if converged
        if(i>0){
            if(abs(this->energy[i-1]-this->energy[i])<this->scf_tol)
                converged=true;
        }
    }

    // end
    F.clear();
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
    orbital_basis_change(basis_change);

    // initialize Helmholtz operator
    S.calculate_helmholtz_coefficients();
    std::shared_ptr<DoubleVector> ptr_helm_coeff = S.get_helmholtz_coefficients();
    HelmholtzVector H(this->prec, *ptr_helm_coeff);
    
    // new orbitals
    OrbitalVector new_Phi(L);

    orbital_update_one_body(F, S, new_Phi);
    orbital_update_two_body(F, S, new_Phi);
    
    // apply Helmholtz kernel and update orbitals
    (*this->orbitals) = H(new_Phi);
}

void LagrangianSolver::orbital_basis_change(std::shared_ptr<ComplexMatrix> basis_change){
    const int L = this->orbitals->size();
    OrbitalVector new_Phi;

    for (int i = 0; i < L; i++) {
        std::vector<ComplexDouble> coeffs(basis_change->row(i).data(), 
                                          basis_change->row(i).data() + L);
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
                    mrcpp::add(new_Phi[m], 1.0, new_Phi[m], -lag_coeff(m, j) / one_rdm(m, m), (*this->orbitals)[j], this->prec);
            }
        } else if (std::abs(lag_coeff(m, m)) > this->threshold) {
            std::cerr << "WARNING: lambda_m below threshold, but not helm_coeff'_mm !" << std::endl;
            for (auto j = 0; j < L; j++) {
                if (j != m)
                    mrcpp::add(new_Phi[m], 1.0, new_Phi[m], -lag_coeff(m, j) / lag_coeff(m, m), (*this->orbitals)[j], this->prec);
            }
        } else {
            std::cerr << "WARNING: Both lambda_m and helm_coeff'_mm below threshold" << std::endl;
            new_Phi[m] = (*this->orbitals)[m];
        }
    }
}

void LagrangianSolver::orbital_update_two_body(FockBuilder &F, ChemTensorSolver &S, OrbitalVector &new_Phi){
    const int L = this->orbitals->size();
    auto &one_rdm = *(S.get_one_rdm());
    auto &two_rdm = *(S.get_two_rdm());
    auto &helm_coeff = *(S.get_helmholtz_coefficients());
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

    // two-body term
    for (int j = 0; j < L; j++) {
        std::cout << "Orbital " << j+1 << " out of " << L << std::endl;
        for (int l = 0; l < L; l++) {
            // g_jl = 4*pi * poisson(Phi[j] * Phi[l])
            // TODO: probably BUG here!
            Orbital phi_jl, g_jl;
            mrcpp::multiply(phi_jl, (*this->orbitals)[j], (*this->orbitals)[l], this->prec);
            mrcpp::apply(this->prec, g_jl, P, phi_jl);

            for (int k = 0; k < L; k++) {
                // state = g_jl * Phi[k]
                Orbital state;
                mrcpp::multiply(state, g_jl, (*this->orbitals)[k], this->prec);

                for (int m = 0; m < L; m++) {
                    ComplexDouble denom;
                    if (std::abs(one_rdm(m, m)) > threshold)
                        denom = one_rdm(m, m);
                    else if (std::abs(helm_coeff(m, m)) > threshold)
                        denom = helm_coeff(m, m);
                    else
                        continue; // both below threshold, skip

                    ComplexDouble coeff = 2.0 * U_dE2(m, j, k, l) / denom;
                    mrcpp::add(new_Phi[m], 1.0, new_Phi[m], coeff, state, this->prec);
                }
            }
        }
    }
}




} // namespace mrchem