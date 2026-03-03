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

// This include must come first due to name clashes

extern "C" {
#include "hamiltonian.h"
#include "qnumber.h"
}

#include "ChemTensorSolver.h"

namespace mrchem {

ChemTensorSolver::ChemTensorSolver(OrbitalVector &Phi, FockBuilder &F, int Ne, int spin) : ExternalSolver() {
    set_integrals(Phi, F);
    set_dense_tensors();
    this->qnum_sector = encode_quantum_number_pair(Ne, spin);
}

void ChemTensorSolver::set_dense_tensors(){
    this->tkin_tensor = new dense_tensor;
    this->tkin_tensor->data  = static_cast<void*>(this->one_body_integrals->data());
    this->tkin_tensor->dim   = new ct_long[2]{this->one_body_integrals->rows(), this->one_body_integrals->cols()};
    this->tkin_tensor->dtype = CT_DOUBLE_COMPLEX;
    this->tkin_tensor->ndim  = 2;

    this->vnuc_tensor = new dense_tensor;
    this->vnuc_tensor->data  = static_cast<void*>(this->two_body_integrals->data());
    this->vnuc_tensor->dim   = new ct_long[4]{this->two_body_integrals->dimension(0), this->two_body_integrals->dimension(1), this->two_body_integrals->dimension(2), this->two_body_integrals->dimension(3)};
    this->vnuc_tensor->dtype = CT_DOUBLE_COMPLEX;
    this->vnuc_tensor->ndim  = 4;
}

void ChemTensorSolver::optimize() {
    if (!this->one_body_integrals || !this->two_body_integrals) {
        MSG_ABORT("Integrals not set. Call set_integrals() before optimize().");
    }

    mpo hamiltonian;
    mpo_assembly assembly;
    construct_molecular_hamiltonian_mpo_assembly(this->tkin_tensor, this->vnuc_tensor, this->optimize_assembly, &assembly);
    mpo_from_assembly(&assembly, &hamiltonian);
    delete_mpo_assembly(&assembly);

    //this->energy = encode_quantum_number_pair(4, 2); // dummy
    //calculate_rdms();
    
}

void ChemTensorSolver::calculate_rdms() {
    *(this->one_rdm) = *(this->one_body_integrals);
    *(this->two_rdm) = *(this->two_body_integrals);
    // dummy (dimensions are the same)
    auto *elements = this->one_body_integrals->data();
}

} // namespace mrchem