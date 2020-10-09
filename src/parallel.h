/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
extern std::vector<int> bankmaster;

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

void allreduce_vector(IntVector &vec, MPI_Comm comm);
void allreduce_vector(DoubleVector &vec, MPI_Comm comm);
void allreduce_vector(ComplexVector &vec, MPI_Comm comm);
void allreduce_matrix(IntMatrix &vec, MPI_Comm comm);
void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm);
void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm);

} // namespace mpi

namespace bank {

struct deposit {
    Orbital *orb;
    double *data; // for pure data arrays
    bool hasdata;
    int datasize;
    int id;     // to identify what is deposited
    int source; // mpi rank from the source of the data
};

struct queue_struct {
    int id;
    std::vector<int> clients;
};

} // namespace bank

class Bank {
public:
    Bank() = default;
    ~Bank();
    void open();
    void close();
    void clear_all(int i = mpi::orb_rank, MPI_Comm comm = mpi::comm_orb);
    void clear(int ix);
    int put_orb(int id, Orbital &orb);
    int get_orb(int id, Orbital &orb, int wait = 0);
    int get_orb_del(int id, Orbital &orb);
    int put_func(int id, QMFunction &func);
    int get_func(int id, QMFunction &func);
    int set_datasize(int datasize);
    int put_data(int id, int size, double *data);
    int get_data(int id, int size, double *data);
    int get_maxtotalsize();

private:
    int const CLOSE_BANK = 1;
    int const CLEAR_BANK = 2;
    int const GET_ORBITAL = 3;
    int const GET_ORBITAL_AND_WAIT = 4;
    int const GET_ORBITAL_AND_DELETE = 5;
    int const SAVE_ORBITAL = 6;
    int const GET_FUNCTION = 7;
    int const SAVE_FUNCTION = 8;
    int const SET_DATASIZE = 9;
    int const GET_DATA = 10;
    int const SAVE_DATA = 11;
    int const GETMAXTOTDATA = 12;
    std::map<int, int> id2ix;
    std::vector<bank::deposit> deposits;
    std::map<int, int> id2qu;
    std::vector<bank::queue_struct> queue;
    long long currentsize = 0; // total deposited data size (without containers)
    long long maxsize = 0;     // max total deposited data size (without containers)

    void clear_bank();
};

namespace mpi {
extern Bank orb_bank;
} // namespace mpi

} // namespace mrchem
