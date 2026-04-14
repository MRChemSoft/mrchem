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
#include "gradient.h"
#include "aligned_memory.h"
}

#include "ChemTensorSolver.h"

namespace mrchem {

ChemTensorSolver::ChemTensorSolver(OrbitalVector &Phi, FockBuilder &F, Nuclei &nucs, int Ne, int spin, json dict_chemtensor) : ExternalSolver(F, nucs) {
    this->qnum_sector = encode_quantum_number_pair(Ne, spin);
    std::cout << "Initialized ChemTensorSolver with Ne = " << Ne << " and spin = " << spin << std::endl;

    this->max_vdim = dict_chemtensor["max_vdim"]; 
    this->num_sweeps = dict_chemtensor["num_sweeps"];
    this->maxiter_lanczos = dict_chemtensor["maxiter_lanczos"];
    this->tol_split = dict_chemtensor["tol_split"];
    this->optimize_assembly = false; //dict_chemtensor["optimize_assembly"];
}

ChemTensorSolver::~ChemTensorSolver() {
    if (this->tkin_tensor) {
        delete this->tkin_tensor->dim;
        delete this->tkin_tensor;
    }
    if (this->velec_tensor) {
        delete this->velec_tensor->dim;
        delete this->velec_tensor;
    }
    if (this->assembly){ 
        delete_mpo_assembly(this->assembly); 
        delete this->assembly; 
    }
    if (this->psi){ 
        delete_mps(this->psi);
        delete this->psi;
    }
}

void ChemTensorSolver::set_dense_tensors(){
    this->tkin_tensor = new dense_tensor;
    this->tkin_tensor->data  = static_cast<void*>(this->one_body_integrals->data());
    this->tkin_tensor->dim   = new ct_long[2]{this->one_body_integrals->rows(), this->one_body_integrals->cols()};
    this->tkin_tensor->dtype = CT_DOUBLE_COMPLEX;
    this->tkin_tensor->ndim  = 2;

    this->velec_tensor = new dense_tensor;
    this->velec_tensor->data  = static_cast<void*>(this->two_body_integrals->data()); //->shuffle(Eigen::array<int,4>{0,2,1,3})); //physicist's notation
    this->velec_tensor->dim   = new ct_long[4]{this->two_body_integrals->dimension(0), this->two_body_integrals->dimension(1), this->two_body_integrals->dimension(2), this->two_body_integrals->dimension(3)};
    this->velec_tensor->dtype = CT_DOUBLE_COMPLEX;
    this->velec_tensor->ndim  = 4;
}

void ChemTensorSolver::set_integrals(OrbitalVector &Phi){
    ExternalSolver::set_integrals(Phi);
    set_dense_tensors();
}

void ChemTensorSolver::optimize() {
    if (!this->one_body_integrals || !this->two_body_integrals)
        MSG_ABORT("Integrals not set.");
    
    if(!this->assembly)
        this->assembly = new mpo_assembly{};
    
    mpo hamiltonian;
    construct_spin_molecular_hamiltonian_mpo_assembly(this->tkin_tensor, this->velec_tensor, this->optimize_assembly, this->assembly);
    //this->assembly = std::make_shared<mpo_assembly>(assembly);

    mpo_from_assembly(this->assembly, &hamiltonian);
    if (!mpo_is_consistent(&hamiltonian))
		MSG_ABORT("internal consistency check for Molecular Hamiltonian MPO failed");
	
    // initial state vector as MPS
    if(!this->psi)
        this->psi = new mps{};
	{
		rng_state rng;
		seed_rng_state(42, &rng);

		construct_random_mps(hamiltonian.a[0].dtype, hamiltonian.nsites, hamiltonian.d, hamiltonian.qsite, this->qnum_sector, this->max_vdim, &rng, this->psi);
		if (!mps_is_consistent(this->psi)) 
			MSG_ABORT("internal MPS consistency check failed");
	}

	// #ifdef _OPENMP
	// printf("maximum number of OpenMP threads: %d\n", omp_get_max_threads());
	// #else
	// printf("OpenMP not available\n");
	// #endif

	// run two-site DMRG
	this->en_sweeps.reserve(this->num_sweeps);
	std::vector<double> entropy(hamiltonian.nsites - 1);
	if (dmrg_twosite(&hamiltonian, this->num_sweeps, this->maxiter_lanczos, this->tol_split, this->max_vdim, this->psi, this->en_sweeps.data(), entropy.data()) < 0)
		std::cerr << "'dmrg_twosite' failed internally" << std::endl;
    
    
    if(this->energy_correction)
        for(auto i=0; i<this->num_sweeps; i++)
            this->en_sweeps.data()[i] += this->E_nn;
    
	this->energy = this->en_sweeps[this->num_sweeps - 1];



    // calculate final bond dimensions
    this->bond_dimensions.reserve(hamiltonian.nsites + 1);
	for (int l = 0; l < hamiltonian.nsites + 1; l++)
		this->bond_dimensions[l] = mps_bond_dim(this->psi, l);
	
    std::cout << "virtual bond dimension: [";
    for (int l = 0; l < hamiltonian.nsites + 1; l++)
        std::cout << this->bond_dimensions[l] << ", ";
    std::cout << "]" << std::endl;

    // calculate RDMs
    calculate_rdms();
}

void ChemTensorSolver::calculate_rdms() {
    // number of orbitals
    const int N = this->psi->nsites;

    // vectors storing energy averages and gradients' coefficients
    std::vector<double> avr(this->assembly->num_coeffs);
    std::vector<std::complex<double>> dcoeff(this->assembly->num_coeffs);
    operator_average_coefficient_gradient(this->assembly, this->psi, this->psi, avr.data(), dcoeff.data());

    // --- one-body RDM ---
    int c = 2;
    ComplexMatrix rdm1(N,N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            rdm1(i, j) = dcoeff[c++];
    // symmetrize
    rdm1 = 0.5 * (rdm1 + rdm1.conjugate().transpose());
    this->one_rdm = std::make_shared<ComplexMatrix>(rdm1);


    // --- two-body RDM ---
    const int dim_g0 = (this->assembly->num_coeffs - c) / 2;

    // intermediate tensors
    ComplexTensorR4 dg0(N, N, N, N); 
    ComplexTensorR4 dg1(N, N, N, N);
    dg0.setZero();
    dg1.setZero();
    int c0 = c;
    int c1 = c + dim_g0;
    for (int i = 0; i < N; i++)
        for (int j = i; j < N; j++)
            for (int k = 0; k < N; k++)
                for (int l = k; l < N; l++) {
                    dg0(i, j, k, l) = dcoeff[c0];
                    dg1(i, j, k, l) = dcoeff[c1];
                    c0++;
                    c1++;
                }

    // reconstruct two_body from dg0 and dg1
    ComplexTensorR4 rdm2(N, N, N, N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                for (int l = 0; l < N; l++)
                    rdm2(i, j, k, l) = 0.5 * (
                          dg0(i, j, k, l) + dg0(j, i, l, k)
                        - dg1(j, i, k, l) - dg1(i, j, l, k));

    // symmetrize: rdm2 + rdm2.conj().transpose(2, 1, 0, 3)
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                for (int l = 0; l < N; l++)
                    rdm2(i, j, k, l) = 0.5 * (rdm2(i, j, k, l) + std::conj(rdm2(k, j, i, l)));
    
    // store two body RDM in variable of parent class
    this->two_rdm = std::make_shared<ComplexTensorR4>(rdm2);
}
} // namespace mrchem