#include "MRCPP/Printer"

#include "parallel.h"
#include "Orbital.h"

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

void initialize(int argc, char **argv) {
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
  
void finalize() {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

} //namespace mpi


/*********************************
 * Orbital related MPI functions *
 *********************************/
namespace orbital {

//send an orbital with MPI
void send_orbital(Orbital &orb, int dst, int tag) {
#ifdef HAVE_MPI
    NEEDS_TESTING;

    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Send(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, dst, 0, mpi::comm_orb);

    if (orb.hasReal()) mrcpp::send_tree(orb.real(), dst, tag, mpi::comm_orb, orbinfo.nChunksReal);
    if (orb.hasImag()) mrcpp::send_tree(orb.imag(), dst, tag+10000, mpi::comm_orb, orbinfo.nChunksImag);
#endif
}

//send an orbital with MPI
void isend_orbital(Orbital &orb, int dst, int tag, MPI_Request& request) {
#ifdef HAVE_MPI
    NEEDS_TESTING;

    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Isend(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, dst, 0, mpi::comm_orb, &request);

    if (orb.hasReal()) mrcpp::isend_tree(orb.real(), dst, tag, mpi::comm_orb,  &request, orbinfo.nChunksReal);
    if (orb.hasImag()) mrcpp::isend_tree(orb.imag(), dst, tag+10000, mpi::comm_orb, &request, orbinfo.nChunksImag);

#endif
}

//receive an orbital with MPI
void recv_orbital(Orbital &orb, int src, int tag) {
#ifdef HAVE_MPI
    NEEDS_TESTING;
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

} //namespace orbital


} //namespace mrchem
