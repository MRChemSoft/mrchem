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

#include "MRCPP/MWOperators"

#include "ExternalSolver.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/MomentumOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/qmoperator_utils.h"
#include "qmoperators/two_electron/two_electron_utils.h"

namespace mrchem {

class FockBuilder;
// class PoissonOperator;

ExternalSolver::ExternalSolver(FockBuilder &F, Nuclei &nucs){
    this->F = F;
    this->nucs = nucs;
    // TODO: read from input
    this->prec = 1e-3;
}


/** @brief Calculates and stores the one- and two-electron integrals for given orbitals
 *
 * @param Phi: Vector of orbitals
 *
 * Calculates the one- and two-electron integrals for the orbitals in Phi, and stores them
 * in the class members one_body_integrals and two_body_integrals.
 *
 */

void ExternalSolver::set_integrals(OrbitalVector &Phi) {
    this->F.setup(this->prec);
    // operators
    MomentumOperator P = this->F.momentum();
    NuclearOperator V = *(this->F.getNuclearOperator());
    
    // set the one- and two-body integrals
    ExternalSolver::set_one_body_integrals(Phi, P, V);
    ExternalSolver::set_two_body_integrals(Phi);

    // set nuclear repulsion energy
    this->E_nn = chemistry::compute_nuclear_repulsion(this->nucs);
}

// Private

// TODO: change 'NuclearOperator' to 'RankZeroOperator'
void ExternalSolver::set_one_body_integrals(OrbitalVector &Phi, MomentumOperator &P, NuclearOperator &V) {
    this->one_body_integrals = std::make_shared<ComplexMatrix>(qmoperator::calc_kinetic_matrix(P, Phi, Phi) + V(Phi, Phi));
}

void ExternalSolver::set_two_body_integrals(OrbitalVector &Phi) {
    this->two_body_integrals = std::make_shared<ComplexTensorR4>(calc_2elintegrals(this->prec, Phi));
}

void ExternalSolver::calculate_lagrange_multipliers(){
    const int L = this->one_body_integrals->rows();
    ComplexMatrix lag_coeff = ComplexMatrix::Zero(L, L);

    // one-body contributions
    lag_coeff += (*this->one_rdm) * (*this->one_body_integrals).transpose();

    // two-body contributions
    for (int n = 0; n < L; n++)
        for (int m = 0; m < L; m++)
            for (int j = 0; j < L; j++)
                for (int k = 0; k < L; k++)
                    for (int l = 0; l < L; l++)
                        lag_coeff(n, m) += 2.0 * (*this->two_rdm)(m, j, k, l) * (*this->two_body_integrals)(n, j, k, l);
    
    this->lag_coeff = std::make_shared<ComplexMatrix>(lag_coeff);
}

void ExternalSolver::diagonalize_1rdm(){
    const int L = this->one_rdm->rows();
    // NOTE: I need to calculate lag_coeff before starting with the basis change!
    if(!this->lag_coeff)
        MSG_ABORT("Lagrange multipliers matrix must be defined before the change of basis.");
        
    if (!this->basis_change)
        this->basis_change = std::make_shared<ComplexMatrix>(ComplexMatrix::Identity(L, L));

    // diagonalize the one-body RDM (Hermitian solver)
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> solver(*this->one_rdm);
    if (solver.info() != Eigen::Success)
        MSG_ABORT("Diagonalization of one-body RDM failed");
    // compute matrix of basis change
    ComplexMatrix eigvecs = solver.eigenvectors();
    // reverse order (largest eigenvalue first)
    ComplexMatrix mat(L, L);
    for (int i = 0; i < L; i++)
        mat.col(i) = eigvecs.col(L - 1 - i);
    // U = mat^* .T so that U @ one_rdm @ U^† = diag(eig_val)
    *(this->basis_change) = mat.adjoint();

    // one-body rdm in new basis
    *this->one_rdm = ComplexMatrix::Zero(L, L);
    for (int i = 0; i < L; i++)
        (*this->one_rdm)(i, i) = solver.eigenvalues()(L - 1 - i);
    // one-body integral in new basis
    *this->one_body_integrals = (*this->basis_change) * (*this->one_body_integrals) * (*this->basis_change).adjoint();
    // lagrange multipliers in the new basis
    *this->lag_coeff = (*this->basis_change) * (*this->lag_coeff) * (*this->basis_change).adjoint();
}

void ExternalSolver::calculate_helmholtz_coefficients() {
    const int L = this->one_body_integrals->rows();

    DoubleVector result(L);
    for (int m = 0; m < L; m++) {
        if (std::abs((*this->one_rdm)(m, m)) < 1e-10)
            MSG_ABORT("Division by zero in helmholtz_coefficients: occupation number too small");
        result[m] = ((*this->lag_coeff)(m, m) / (*this->one_rdm)(m, m)).real();
    }

    this->helm_coeff = std::make_shared<DoubleVector>(result);
}

// std::vector<ComplexDouble> LagrangianSolver::calculate_helmholtz_coefficients() {
//     const int L = this->one_body_integrals->rows();

//     // 2 * np.einsum("mn, mi, ijkl, njkl -> m", U, U, two_rdm, two_body_int)
//     std::vector<ComplexDouble> correction(L, 0.0);
//     for (int m = 0; m < L; m++)
//         for (int n = 0; n < L; n++)
//             for (int i = 0; i < L; i++)
//                 for (int j = 0; j < L; j++)
//                     for (int k = 0; k < L; k++)
//                         for (int l = 0; l < L; l++)
//                             correction[m] += 2.0 
//                                 * (*this->change_basis)(m, n) 
//                                 * (*this->change_basis)(m, i)
//                                 * (*this->two_rdm)(i, j, k, l)
//                                 * (*this->two_body_integrals)(n, j, k, l);

//     std::vector<ComplexDouble> result(L);
//     for (int m = 0; m < L; m++) {
//         if (std::abs((*this->one_rdm)(m, m)) < 1e-10)
//             MSG_ABORT("Division by zero in helmholtz_coefficients: occupation number too small");
//         result[m] = (*this->one_body_integrals)(m, m) + correction[m] / (*this->one_rdm)(m, m);
//     }
//     this->helm_coeff = std::make_shared<ComplexMatrix>(result);
// }



} // namespace mrchem