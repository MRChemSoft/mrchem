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

void ExternalSolver::diagonalize_1rdm(){
    const int L = this->one_rdm->rows();

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

}




} // namespace mrchem