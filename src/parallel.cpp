#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "Orbital.h"
#include "orbital_utils.h"
#include "Density.h"
#include "density_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {


namespace omp {

int n_threads = omp_get_max_threads();

} //namespace omp


namespace mpi {

int orb_rank = 0;
int orb_size = 1;
int share_rank = 0;
int share_size = 1;

MPI_Comm comm_orb;
MPI_Comm comm_share;

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
    MPI_Comm comm_sh_group;
    MPI_Comm_split(MPI_COMM_WORLD, share_rank, world_rank, &comm_sh_group);
    // mpiShRank is color (same color->in same group)
    // MPI_worldrank is key (orders rank within the groups)

    // we define a new orbital rank, so that the orbitals within
    // a shared memory group, have consecutive ranks
    int sh_group_rank;
    MPI_Comm_rank(comm_sh_group, &sh_group_rank);

    orb_rank = share_rank + sh_group_rank*world_size;
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
        if (mpi::my_orb(Phi[i])) {
            chunk.push_back(std::make_tuple(i, Phi[i]));
        }
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

//send an orbital with MPI
void mpi::send_orbital(Orbital &orb, int dst, int tag) {
#ifdef HAVE_MPI
    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Send(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, dst, 0, mpi::comm_orb);

    if (orb.hasReal()) mrcpp::send_tree(orb.real(), dst, tag, mpi::comm_orb, orbinfo.nChunksReal);
    if (orb.hasImag()) mrcpp::send_tree(orb.imag(), dst, tag+10000, mpi::comm_orb, orbinfo.nChunksImag);
#endif
}

//send an orbital with MPI
void mpi::isend_orbital(Orbital &orb, int dst, int tag, MPI_Request& request) {
#ifdef HAVE_MPI
    NEEDS_TESTING;

    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Isend(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, dst, 0, mpi::comm_orb, &request);

    if (orb.hasReal()) mrcpp::isend_tree(orb.real(), dst, tag, mpi::comm_orb,  &request, orbinfo.nChunksReal);
    if (orb.hasImag()) mrcpp::isend_tree(orb.imag(), dst, tag+10000, mpi::comm_orb, &request, orbinfo.nChunksImag);

#endif
}

//receive an orbital with MPI
void mpi::recv_orbital(Orbital &orb, int src, int tag) {
#ifdef HAVE_MPI
    MPI_Status status;

    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Recv(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, src, 0, mpi::comm_orb, &status);

    if (orbinfo.nChunksReal > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (orb.hasReal()) MSG_FATAL("Orbital not empty");
        orb.alloc(NUMBER::Real);
        mrcpp::recv_tree(orb.real(), src, tag, mpi::comm_orb, orbinfo.nChunksReal);
    }

    if (orbinfo.nChunksImag > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (orb.hasImag()) MSG_FATAL("Orbital not empty");
        orb.alloc(NUMBER::Imag);
        mrcpp::recv_tree(orb.imag(), src, tag+10000, mpi::comm_orb, orbinfo.nChunksImag);
    }
#endif
}

void mpi::reduce_density(Density &rho, MPI_Comm comm) {
#ifdef HAVE_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    if (comm_size == 1) return;

    if (not rho.hasReal()) NOT_IMPLEMENTED_ABORT;
    if (rho.hasImag()) NOT_IMPLEMENTED_ABORT;

    Timer timer;
    if (comm_rank == 0) {
        Density *rho_0 = new Density;
        rho_0->alloc(NUMBER::Real);
        mrcpp::copy_grid(rho_0->real(), rho.real());
        mrcpp::copy_func(rho_0->real(), rho.real());
        rho.real().clear();

        mrcpp::FunctionTreeVector<3> rho_vec;
        rho_vec.push_back(std::make_tuple(1.0, &rho_0->real()));
        for (int src = 1; src < comm_size; src++) {
            Density *rho_i = new Density;
            int tag = 3333+src;
            rho_i->alloc(NUMBER::Real);
            mrcpp::recv_tree(rho_i->real(), src, tag, comm);
            if (rho_i->real().getSquareNorm() > 0.0) rho_vec.push_back(std::make_tuple(1.0, &rho_i->real()));
        }
        mrcpp::build_grid(rho.real(), rho_vec);
        mrcpp::add(-1.0, rho.real(), rho_vec, 0);
        mrcpp::clear(rho_vec, true);
    } else {
        int tag = 3333+comm_rank;
        mrcpp::send_tree(rho.real(), 0, tag, comm);
        rho.real().clear();
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

    if (not rho.hasReal()) NOT_IMPLEMENTED_ABORT;
    if (rho.hasImag()) NOT_IMPLEMENTED_ABORT;

    Timer timer;
    if (comm_rank == 0) {
        for (int dst = 1; dst < comm_size; dst++) {
            int tag = 4444+dst;
            mrcpp::send_tree(rho.real(), dst, tag, comm);
        }
    } else {
        int tag = 4444+comm_rank;
        mrcpp::recv_tree(rho.real(), 0, tag, comm);
    }
    MPI_Barrier(comm);
    timer.stop();
    Printer::printDouble(0, "Broadcast density", timer.getWallTime(), 5);
#endif
}

} //namespace mrchem
