#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/trees/FunctionNode.h"

#include "parallel.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace omp {

int n_threads = omp_get_max_threads();

} //namespace omp

namespace mpi {

bool numerically_exact = false;
bool share_nuc_pot = false;
bool share_coul_dens = false;
bool share_coul_pot = false;
bool share_xc_dens = false;
bool share_xc_pot = false;
int shared_memory_size = 0;

int orb_rank = 0;
int orb_size = 1;
int share_rank = 0;
int share_size = 1;
int sh_group_rank = 0;

MPI_Comm comm_orb;
MPI_Comm comm_share;
MPI_Comm comm_sh_group;

} //namespace mpi

void mpi::initialize(int argc, char **argv) {
    omp_set_dynamic(0);

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // divide the world into groups
    // each group has its own group communicator definition

    // split world into groups that can share memory
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_share);

    MPI_Comm_rank(comm_share, &share_rank);
    MPI_Comm_size(comm_share, &share_size);

    // define a rank of the group
    MPI_Comm_split(MPI_COMM_WORLD, share_rank, world_rank, &comm_sh_group);
    // mpiShRank is color (same color->in same group)
    // MPI_worldrank is key (orders rank within the groups)

    // we define a new orbital rank, so that the orbitals within
    // a shared memory group, have consecutive ranks
    MPI_Comm_rank(comm_sh_group, &sh_group_rank);

    orb_rank = share_rank + sh_group_rank * world_size;
    MPI_Comm_split(MPI_COMM_WORLD, 0, orb_rank, &comm_orb);
    // 0 is color (same color->in same group)
    // mpiOrbRank is key (orders rank in the group)

    MPI_Comm_rank(comm_orb, &orb_rank);
    MPI_Comm_size(comm_orb, &orb_size);
#endif
}

void mpi::finalize() {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

void mpi::barrier(MPI_Comm comm) {
#ifdef HAVE_MPI
    MPI_Barrier(comm);
#endif
}

/*********************************
 * Orbital related MPI functions *
 *********************************/

bool mpi::grand_master() {
    return (mpi::orb_rank == 0) ? true : false;
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

/** @brief Clear all orbitals not belonging to this MPI rank */
void mpi::free_foreign(OrbitalVector &Phi) {
    for (int i = 0; i < Phi.size(); i++) {
        if (not mpi::my_orb(Phi[i])) Phi[i].free();
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
#ifdef HAVE_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(DoubleVector &vec, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(ComplexVector &vec, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(IntMatrix &mat, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

//send an orbital with MPI, includes orbital meta data
void mpi::send_orbital(Orbital &orb, int dst, int tag) {
#ifdef HAVE_MPI
    mpi::send_function(orb, dst, tag, mpi::comm_orb);

    OrbitalData &orbinfo = orb.getOrbitalData();
    MPI_Send(&orbinfo, sizeof(OrbitalData), MPI_BYTE, dst, 0, mpi::comm_orb);
#endif
}

//receive an orbital with MPI, includes orbital meta data
void mpi::recv_orbital(Orbital &orb, int src, int tag) {
#ifdef HAVE_MPI
    mpi::recv_function(orb, src, tag, mpi::comm_orb);

    MPI_Status status;
    OrbitalData &orbinfo = orb.getOrbitalData();
    MPI_Recv(&orbinfo, sizeof(OrbitalData), MPI_BYTE, src, 0, mpi::comm_orb, &status);
#endif
}

//send a density with MPI
void mpi::send_function(QMFunction &func, int dst, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    FunctionData &funcinfo = func.getFunctionData();
    MPI_Send(&funcinfo, sizeof(FunctionData), MPI_BYTE, dst, 0, comm);

    if (func.hasReal()) mrcpp::send_tree(func.real(), dst, tag, comm, funcinfo.nChunksReal);
    if (func.hasImag()) mrcpp::send_tree(func.imag(), dst, tag + 10000, comm, funcinfo.nChunksImag);
#endif
}

//receive a denstity with MPI
void mpi::recv_function(QMFunction &func, int src, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    MPI_Status status;

    FunctionData &funcinfo = func.getFunctionData();
    MPI_Recv(&funcinfo, sizeof(FunctionData), MPI_BYTE, src, 0, comm, &status);

    if (funcinfo.nChunksReal > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasReal()) func.alloc(NUMBER::Real);
        mrcpp::recv_tree(func.real(), src, tag, comm, funcinfo.nChunksReal);
    }

    if (funcinfo.nChunksImag > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasImag()) func.alloc(NUMBER::Imag);
        mrcpp::recv_tree(func.imag(), src, tag + 10000, comm, funcinfo.nChunksImag);
    }
#endif
}

/** Update a shared function after it has been changed by one of the MPI ranks. */
void mpi::share_function(QMFunction &func, int src, int tag) {
#ifdef HAVE_MPI
    if (func.isShared()) {
        if (func.hasReal()) mrcpp::share_tree(func.real(), src, tag, mpi::comm_share);
        if (func.hasImag()) mrcpp::share_tree(func.imag(), src, 2 * tag, mpi::comm_share);
    }
#endif
}

/** @brief Add all mpi densities in rank zero density */
void mpi::reduce_density(double prec, Density &rho, MPI_Comm comm) {
#ifdef HAVE_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    if (comm_size == 1) return;

    Timer timer;
    if (comm_rank == 0) {
        for (int src = 1; src < comm_size; src++) {
            Density rho_i(false);
            int tag = 3333 + src;
            mpi::recv_function(rho_i, src, tag, comm);
            rho.add(1.0, rho_i); // add in place using union grid
            rho.crop(prec);
            rho_i.free();
        }
    } else {
        int tag = 3333 + comm_rank;
        mpi::send_function(rho, 0, tag, comm);
    }
    MPI_Barrier(comm);
    timer.stop();
    Printer::printDouble(0, "Reduce density", timer.getWallTime(), 5);
#endif
}

void mpi::broadcast_density(Density &rho, MPI_Comm comm) {
#ifdef HAVE_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    if (comm_size == 1) return;

    Timer timer;
    if (comm_rank == 0) {
        for (int dst = 1; dst < comm_size; dst++) {
            int tag = 4334 + dst;
            mpi::send_function(rho, dst, tag, comm);
        }
    } else {
        int tag = 4334 + comm_rank;
        mpi::recv_function(rho, 0, tag, comm);
    }
    MPI_Barrier(comm);
    timer.stop();
    Printer::printDouble(0, "Broadcast density", timer.getWallTime(), 5);
#endif
}

} //namespace mrchem
