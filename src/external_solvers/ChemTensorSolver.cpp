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
#include "mps.h"
#include "rng.h"
#include "dmrg.h"
#include "qnumber.h"
}

#include "ChemTensorSolver.h"

namespace mrchem {

ChemTensorSolver::ChemTensorSolver(OrbitalVector &Phi, FockBuilder &F, int Ne, int spin, json dict_chemtensor) : ExternalSolver() {
    set_integrals(Phi, F);
    set_dense_tensors();
    this->qnum_sector = encode_quantum_number_pair(Ne, spin);
    std::cout << "Initialized ChemTensorSolver with Ne = " << Ne << " and spin = " << spin << std::endl;

    this->max_vdim = dict_chemtensor["max_vdim"]; 
    this->num_sweeps = dict_chemtensor["num_sweeps"];
    this->maxiter_lanczos = dict_chemtensor["maxiter_lanczos"];
    this->tol_split = dict_chemtensor["tol_split"];
    this->optimize_assembly = dict_chemtensor["optimize_assembly"];

    this->bond_dimensions = nullptr;
    this->en_sweeps = nullptr;
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
        MSG_ABORT("Integrals not set.");
    }

    // overall quantum number sector of quantum state (particle number and spin)

    mpo hamiltonian;
    {
        mpo_assembly assembly;
        construct_molecular_hamiltonian_mpo_assembly(this->tkin_tensor, this->vnuc_tensor, this->optimize_assembly, &assembly);
        mpo_from_assembly(&assembly, &hamiltonian);
        delete_mpo_assembly(&assembly);
    }
    if (!mpo_is_consistent(&hamiltonian))
		MSG_ABORT("internal consistency check for Molecular Hamiltonian MPO failed");
	

    // initial state vector as MPS
	mps psi;
	{
		rng_state rng;
		seed_rng_state(42, &rng);

		construct_random_mps(hamiltonian.a[0].dtype, hamiltonian.nsites, hamiltonian.d, hamiltonian.qsite, this->qnum_sector, this->max_vdim, &rng, &psi);
		if (!mps_is_consistent(&psi)) 
			MSG_ABORT("internal MPS consistency check failed");
		
	}

	// #ifdef _OPENMP
	// printf("maximum number of OpenMP threads: %d\n", omp_get_max_threads());
	// #else
	// printf("OpenMP not available\n");
	// #endif

	// run two-site DMRG
	this->en_sweeps = new double[this->num_sweeps];
	double* entropy   = new double[hamiltonian.nsites - 1];
	if (dmrg_twosite(&hamiltonian, this->num_sweeps, this->maxiter_lanczos, this->tol_split, this->max_vdim, &psi, this->en_sweeps, entropy) < 0)
		std::cerr << "'dmrg_twosite' failed internally" << std::endl;
	this->energy = this->en_sweeps[this->num_sweeps - 1];
    delete [] entropy;

    // calculate final bond dimensions
    this->bond_dimensions = new int[hamiltonian.nsites + 1];
	for (int l = 0; l < hamiltonian.nsites + 1; l++) {
		this->bond_dimensions[l] = mps_bond_dim(&psi, l);
	}





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