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

#pragma once

#include "MRCPP/Parallel"
#include <map>

#ifdef MRCHEM_HAS_MPI
#ifndef MRCPP_HAS_MPI
#include <mpi.h>
#endif
#else
#ifndef MRCPP_HAS_MPI
using MPI_Comm = int;
#endif
#endif

#include "mrchem.h"
#include "qmfunctions/qmfunction_fwd.h"

namespace mrchem {

namespace omp {
extern int n_threads;
} // namespace omp

class Bank;
extern Bank dataBank;

namespace mpi {

extern bool numerically_exact;
extern int shared_memory_size;

extern int world_rank;
extern int world_size;
extern int orb_rank;
extern int orb_size;
extern int share_rank;
extern int share_size;
extern int sh_group_rank;
extern int is_bank;
extern int is_bankclient;
extern int bank_size;
extern int tot_bank_size;
extern int max_tag;
extern std::vector<int> bankmaster;
extern int task_bank;

extern MPI_Comm comm_orb;
extern MPI_Comm comm_share;
extern MPI_Comm comm_sh_group;
extern MPI_Comm comm_bank;

void initialize();
void finalize();
void barrier(MPI_Comm comm);

bool grand_master();
bool share_master();
bool my_orb(const Orbital &orb);
bool my_unique_orb(const Orbital &orb);
void distribute(OrbitalVector &Phi);
void free_foreign(OrbitalVector &Phi);
OrbitalChunk get_my_chunk(OrbitalVector &Phi);

void send_orbital(Orbital &orb, int dst, int tag, MPI_Comm comm = mpi::comm_orb);
void recv_orbital(Orbital &orb, int src, int tag, MPI_Comm comm = mpi::comm_orb);

void send_function(QMFunction &func, int dst, int tag, MPI_Comm comm);
void recv_function(QMFunction &func, int src, int tag, MPI_Comm comm);
void share_function(QMFunction &func, int src, int tag, MPI_Comm comm);

void reduce_function(double prec, QMFunction &func, MPI_Comm comm);
void broadcast_function(QMFunction &func, MPI_Comm comm);

void reduce_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, MPI_Comm comm);
void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, OrbitalVector &Phi, MPI_Comm comm);
void broadcast_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, MPI_Comm comm);

void allreduce_vector(IntVector &vec, MPI_Comm comm);
void allreduce_vector(DoubleVector &vec, MPI_Comm comm);
void allreduce_vector(ComplexVector &vec, MPI_Comm comm);
void allreduce_matrix(IntMatrix &vec, MPI_Comm comm);
void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm);
void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm);

} // namespace mpi

} // namespace mrchem
