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

#include "parallel.h"
#include "qmfunctions/ComplexFunction.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "utils/Bank.h"

#ifdef MRCHEM_HAS_OMP
#ifndef MRCPP_HAS_OMP
#include <omp.h>
#endif
#define mrchem_get_max_threads() omp_get_max_threads()
#define mrchem_get_num_threads() omp_get_num_threads()
#define mrchem_get_thread_num() omp_get_thread_num()
#define mrchem_set_dynamic(n) omp_set_dynamic(n)
#else
#define mrchem_get_max_threads() 1
#define mrchem_get_num_threads() 1
#define mrchem_get_thread_num() 0
#define mrchem_set_dynamic(n)
#endif

using mrcpp::Printer;

namespace mrchem {

namespace omp {

int n_threads = mrchem_get_max_threads();

} // namespace omp

using namespace Eigen;

Bank dataBank;

namespace mpi {

bool numerically_exact = false;
int shared_memory_size = 1000;

// these parameters set by initialize()
int world_size = 1;
int world_rank = 0;
int orb_size = 1;
int orb_rank = 0;
int share_size = 1;
int share_rank = 0;
int sh_group_rank = 0;
int is_bank = 0;
int is_centralbank = 0;
int is_bankclient = 1;
int is_bankmaster = 0; // only one bankmaster is_bankmaster
int bank_size = 0;
int tot_bank_size = 0; // size of bank, including the task manager
int max_tag = 0;       // max value allowed by MPI
std::vector<int> bankmaster;
int task_bank = -1; // world rank of the task manager

MPI_Comm comm_orb;
MPI_Comm comm_share;
MPI_Comm comm_sh_group;
MPI_Comm comm_bank;

} // namespace mpi

int id_shift; // to ensure that nodes, orbitals and functions do not collide

extern int metadata_block[3]; // can add more metadata in future
extern int const size_metadata = 3;

void mpi::initialize() {
    Eigen::setNbThreads(1);
    mrchem_set_dynamic(0);
    mrcpp::set_max_threads(omp::n_threads);

#ifdef MRCHEM_HAS_MPI
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::world_rank);

    // divide the world into groups
    // each group has its own group communicator definition

    // define independent group of MPI processes, that are not part of comm_orb
    // for now the new group does not include comm_share
    mpi::comm_bank = MPI_COMM_WORLD; // clients and master
    MPI_Comm comm_remainder;         // clients only

    // set bank_size automatically if not defined by user
    if (mpi::world_size < 2) {
        mpi::bank_size = 0;
    } else if (mpi::bank_size < 0) {
        mpi::bank_size = std::max(mpi::world_size / 3, 1);
    }
    if (mpi::world_size - mpi::bank_size < 1) MSG_ABORT("No MPI ranks left for working!");
    if (mpi::bank_size < 1 and mpi::world_size > 1) MSG_ABORT("Bank size must be at least one when using MPI!");

    mpi::bankmaster.resize(mpi::bank_size);
    for (int i = 0; i < mpi::bank_size; i++) {
        mpi::bankmaster[i] = mpi::world_size - i - 1; // rank of the bankmasters
    }
    if (mpi::world_rank < mpi::world_size - mpi::bank_size) {
        // everything which is left
        mpi::is_bank = 0;
        mpi::is_centralbank = 0;
        mpi::is_bankclient = 1;
    } else {
        // special group of centralbankmasters
        mpi::is_bank = 1;
        mpi::is_centralbank = 1;
        mpi::is_bankclient = 0;
        if (mpi::world_rank == mpi::world_size - mpi::bank_size) mpi::is_bankmaster = 1;
    }
    MPI_Comm_split(MPI_COMM_WORLD, mpi::is_bankclient, mpi::world_rank, &comm_remainder);

    // split world into groups that can share memory
    MPI_Comm_split_type(comm_remainder, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &mpi::comm_share);

    MPI_Comm_rank(mpi::comm_share, &mpi::share_rank);
    MPI_Comm_size(mpi::comm_share, &mpi::share_size);

    // define a rank of the group
    MPI_Comm_split(comm_remainder, mpi::share_rank, mpi::world_rank, &mpi::comm_sh_group);
    // mpiShRank is color (same color->in same group)
    // MPI_worldrank is key (orders rank within the groups)

    // we define a new orbital rank, so that the orbitals within
    // a shared memory group, have consecutive ranks
    MPI_Comm_rank(mpi::comm_sh_group, &mpi::sh_group_rank);

    mpi::orb_rank = mpi::share_rank + mpi::sh_group_rank * mpi::world_size;
    MPI_Comm_split(comm_remainder, 0, mpi::orb_rank, &mpi::comm_orb);
    // 0 is color (same color->in same group)
    // mpiOrbRank is key (orders rank in the group)

    MPI_Comm_rank(mpi::comm_orb, &mpi::orb_rank);
    MPI_Comm_size(mpi::comm_orb, &mpi::orb_size);

    // if bank_size is large enough, we reserve one as "task manager"
    mpi::tot_bank_size = mpi::bank_size;
    if (mpi::bank_size <= 2 and mpi::bank_size > 0) {
        // use the first bank as task manager
        mpi::task_bank = mpi::bankmaster[0];
    } else if (mpi::bank_size > 1) {
        // reserve one bank for task management only
        mpi::bank_size--;
        mpi::task_bank = mpi::bankmaster[mpi::bank_size]; // the last rank is reserved as task manager
    }

    // determine the maximum value alowed for mpi tags
    void *val;
    int flag;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &val, &flag); // max value allowed by MPI for tags
    max_tag = *(int *)val / 2;
    id_shift = max_tag / 2; // half is reserved for non orbital.

    if (mpi::is_bank) {
        // bank is open until end of program
        if (mpi::is_centralbank) { dataBank.open(); }
        mpi::finalize();
        exit(EXIT_SUCCESS);
    }
#else
    mpi::bank_size = 0;
#endif
}

void mpi::finalize() {
#ifdef MRCHEM_HAS_MPI
    if (mpi::bank_size > 0 and mpi::grand_master()) {
        println(4, " max data in bank " << dataBank.get_maxtotalsize() << " MB ");
        dataBank.close();
    }
    MPI_Barrier(MPI_COMM_WORLD); // to ensure everybody got here
    MPI_Finalize();
#endif
}

void mpi::barrier(MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    MPI_Barrier(comm);
#endif
}

/*********************************
 * Orbital related MPI functions *
 *********************************/

bool mpi::grand_master() {
    return (mpi::world_rank == 0 and is_bankclient) ? true : false;
}

bool mpi::share_master() {
    return (mpi::share_rank == 0) ? true : false;
}

/** @brief Test if orbital belongs to this MPI rank (or is common)*/
bool mpi::my_orb(const Orbital &orb) {
    return (orb.rankID() < 0 or orb.rankID() == mpi::orb_rank) ? true : false;
}

/** @brief Test if orbital belongs to this MPI rank */
bool mpi::my_unique_orb(const Orbital &orb) {
    return (orb.rankID() == mpi::orb_rank) ? true : false;
}

/** @brief Distribute orbitals in vector round robin. Orbitals should be empty.*/
void mpi::distribute(OrbitalVector &Phi) {
    for (int i = 0; i < Phi.size(); i++) Phi[i].setRankID(i % mpi::orb_size);
}

/** @brief Free all function pointers not belonging to this MPI rank */
void mpi::free_foreign(OrbitalVector &Phi) {
    for (auto &i : Phi) {
        if (not mpi::my_orb(i)) i.free(NUMBER::Total);
    }
}

/** @brief Return the subset of an OrbitalVector that belongs to this MPI rank */
OrbitalChunk mpi::get_my_chunk(OrbitalVector &Phi) {
    OrbitalChunk chunk;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) chunk.push_back(std::make_tuple(i, Phi[i]));
    }
    return chunk;
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(IntVector &vec, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(DoubleVector &vec, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(ComplexVector &vec, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(IntMatrix &mat, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

// send an orbital with MPI, includes orbital meta data
void mpi::send_orbital(Orbital &orb, int dst, int tag, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    mpi::send_function(orb, dst, tag, comm);
    OrbitalData &orbinfo = orb.getOrbitalData();
    MPI_Send(&orbinfo, sizeof(OrbitalData), MPI_BYTE, dst, 0, comm);
#endif
}

// receive an orbital with MPI, includes orbital meta data
void mpi::recv_orbital(Orbital &orb, int src, int tag, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    mpi::recv_function(orb, src, tag, comm);

    MPI_Status status;
    OrbitalData &orbinfo = orb.getOrbitalData();
    MPI_Recv(&orbinfo, sizeof(OrbitalData), MPI_BYTE, src, 0, comm, &status);
#endif
}

// send a function with MPI
void mpi::send_function(QMFunction &func, int dst, int tag, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
#ifdef MRCPP_HAS_MPI
    if (func.isShared()) MSG_WARN("Sending a shared function is not recommended");
    FunctionData &funcinfo = func.getFunctionData();
    MPI_Send(&funcinfo, sizeof(FunctionData), MPI_BYTE, dst, 0, comm);
    if (func.hasReal()) mrcpp::send_tree(func.real(), dst, tag, comm, funcinfo.real_size);
    if (func.hasImag()) mrcpp::send_tree(func.imag(), dst, tag + 10000, comm, funcinfo.imag_size);
#else
    MSG_ABORT("MRCPP compiled without MPI support");
#endif
#endif
}

// receive a function with MPI
void mpi::recv_function(QMFunction &func, int src, int tag, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
#ifdef MRCPP_HAS_MPI
    if (func.isShared()) MSG_WARN("Receiving a shared function is not recommended");
    MPI_Status status;

    FunctionData &funcinfo = func.getFunctionData();
    MPI_Recv(&funcinfo, sizeof(FunctionData), MPI_BYTE, src, 0, comm, &status);
    if (funcinfo.real_size > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasReal()) func.alloc(NUMBER::Real);
        mrcpp::recv_tree(func.real(), src, tag, comm, funcinfo.real_size);
    }

    if (funcinfo.imag_size > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasImag()) func.alloc(NUMBER::Imag);
        mrcpp::recv_tree(func.imag(), src, tag + 10000, comm, funcinfo.imag_size);
    }
#else
    MSG_ABORT("MRCPP compiled without MPI support");
#endif
#endif
}

/** Update a shared function after it has been changed by one of the MPI ranks. */
void mpi::share_function(QMFunction &func, int src, int tag, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
#ifdef MRCPP_HAS_MPI
    if (func.isShared()) {
        if (func.hasReal()) mrcpp::share_tree(func.real(), src, tag, comm);
        if (func.hasImag()) mrcpp::share_tree(func.imag(), src, 2 * tag, comm);
    }
#else
    MSG_ABORT("MRCPP compiled without MPI support");
#endif
#endif
}

/** @brief Add all mpi function into rank zero */
void mpi::reduce_function(double prec, QMFunction &func, MPI_Comm comm) {
/* 1) Each odd rank send to the left rank
   2) All odd ranks are "deleted" (can exit routine)
   3) new "effective" ranks are defined within the non-deleted ranks
      effective rank = rank/fac , where fac are powers of 2
   4) repeat
 */
#ifdef MRCHEM_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) {
        if ((comm_rank / fac) % 2 == 0) {
            // receive
            int src = comm_rank + fac;
            if (src < comm_size) {
                QMFunction func_i(false);
                int tag = 3333 + src;
                mpi::recv_function(func_i, src, tag, comm);
                func.add(1.0, func_i); // add in place using union grid
                func.crop(prec);
            }
        }
        if ((comm_rank / fac) % 2 == 1) {
            // send
            int dest = comm_rank - fac;
            if (dest >= 0) {
                int tag = 3333 + comm_rank;
                mpi::send_function(func, dest, tag, comm);
                break; // once data is sent we are done
            }
        }
        fac *= 2;
    }
    MPI_Barrier(comm);
#endif
}

/** @brief make union tree and send into rank zero */
void mpi::reduce_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, MPI_Comm comm) {
/* 1) Each odd rank send to the left rank
   2) All odd ranks are "deleted" (can exit routine)
   3) new "effective" ranks are defined within the non-deleted ranks
      effective rank = rank/fac , where fac are powers of 2
   4) repeat
 */
#ifdef MRCHEM_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) {
        if ((comm_rank / fac) % 2 == 0) {
            // receive
            int src = comm_rank + fac;
            if (src < comm_size) {
                int tag = 3333 + src;
                mrcpp::FunctionTree<3> tree_i(*MRA);
                mrcpp::recv_tree(tree_i, src, tag, comm, -1, false);
                tree.appendTreeNoCoeff(tree_i); // make union grid
            }
        }
        if ((comm_rank / fac) % 2 == 1) {
            // send
            int dest = comm_rank - fac;
            if (dest >= 0) {
                int tag = 3333 + comm_rank;
                mrcpp::send_tree(tree, dest, tag, comm, -1, false);
                break; // once data is sent we are done
            }
        }
        fac *= 2;
    }
    MPI_Barrier(comm);
#endif
}

/** @brief make union tree without coeff and send to all
 *  Include both real and imaginary parts
 */
void mpi::allreduce_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, OrbitalVector &Phi, MPI_Comm comm) {
    /* 1) make union grid of own orbitals
       2) make union grid with others orbitals (sent to rank zero)
       3) rank zero broadcast func to everybody
     */
    int N = Phi.size();
    for (int j = 0; j < N; j++) {
        if (not mpi::my_orb(Phi[j])) continue;
        if (Phi[j].hasReal()) tree.appendTreeNoCoeff(Phi[j].real());
        if (Phi[j].hasImag()) tree.appendTreeNoCoeff(Phi[j].imag());
    }
#ifdef MRCHEM_HAS_MPI
    mpi::reduce_Tree_noCoeff(tree, mpi::comm_orb);
    mpi::broadcast_Tree_noCoeff(tree, mpi::comm_orb);
#endif
}

/** @brief Distribute rank zero function to all ranks */
void mpi::broadcast_function(QMFunction &func, MPI_Comm comm) {
/* use same strategy as a reduce, but in reverse order */
#ifdef MRCHEM_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) fac *= 2;
    fac /= 2;

    while (fac > 0) {
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 1) {
            // receive
            int src = comm_rank - fac;
            int tag = 4334 + comm_rank;
            mpi::recv_function(func, src, tag, comm);
        }
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 0) {
            // send
            int dst = comm_rank + fac;
            int tag = 4334 + dst;
            if (dst < comm_size) mpi::send_function(func, dst, tag, comm);
        }
        fac /= 2;
    }
    MPI_Barrier(comm);
#endif
}

/** @brief Distribute rank zero function to all ranks */
void mpi::broadcast_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, MPI_Comm comm) {
/* use same strategy as a reduce, but in reverse order */
#ifdef MRCHEM_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) fac *= 2;
    fac /= 2;

    while (fac > 0) {
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 1) {
            // receive
            int src = comm_rank - fac;
            int tag = 4334 + comm_rank;
            mrcpp::recv_tree(tree, src, tag, comm, -1, false);
        }
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 0) {
            // send
            int dst = comm_rank + fac;
            int tag = 4334 + dst;
            if (dst < comm_size) mrcpp::send_tree(tree, dst, tag, comm, -1, false);
        }
        fac /= 2;
    }
    MPI_Barrier(comm);
#endif
}

} // namespace mrchem
