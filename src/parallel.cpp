#include "MRCPP/Printer"

#include "parallel.h"
#include "Orbital.h"

namespace mrchem {

int mpiOrbRank = 0;
int mpiOrbSize = 1;
int mpiShRank = 0;
int mpiShSize = 1;

MPI_Comm mpiCommOrb;
MPI_Comm mpiCommSh;

namespace orbital {
/*********************************
 * Orbital related MPI functions *
 *********************************/

//send an orbital with MPI
void send_orbital(Orbital &orb, int dst, int tag) {
#ifdef HAVE_MPI
    NEEDS_TESTING;
    MPI_Comm comm = mpiCommOrb;

    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Send(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, dst, 0, comm);

    if (orb.hasReal()) mrcpp::send_tree(orb.real(), dst, tag, comm, orbinfo.nChunksReal);
    if (orb.hasImag()) mrcpp::send_tree(orb.imag(), dst, tag+10000, comm, orbinfo.nChunksImag);
#endif
}

//send an orbital with MPI
void isend_orbital(Orbital &orb, int dst, int tag, MPI_Request& request) {
#ifdef HAVE_MPI
    NEEDS_TESTING;
    MPI_Comm comm = mpiCommOrb;

    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Isend(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, dst, 0, comm, &request);

    if (orb.hasReal()) mrcpp::isend_tree(orb.real(), dst, tag, comm,  &request, orbinfo.nChunksReal);
    if (orb.hasImag()) mrcpp::isend_tree(orb.imag(), dst, tag+10000, comm, &request, orbinfo.nChunksImag);

#endif
}

//receive an orbital with MPI
void recv_orbital(Orbital &orb, int src, int tag) {
#ifdef HAVE_MPI
    NEEDS_TESTING;
    MPI_Status status;
    MPI_Comm comm = mpiCommOrb;

    OrbitalMeta &orbinfo = orb.getMetaData();
    MPI_Recv(&orbinfo, sizeof(OrbitalMeta), MPI_BYTE, src, 0, comm, &status);

    if (orbinfo.nChunksReal > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (orb.hasReal()) MSG_FATAL("Orbital not empty");
        orb.alloc(NUMBER::Real);
        mrcpp::recv_tree(orb.real(), src, tag, comm, orbinfo.nChunksReal);
    }

    if (orbinfo.nChunksImag > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (orb.hasImag()) MSG_FATAL("Orbital not empty");
        orb.alloc(NUMBER::Imag);
        mrcpp::recv_tree(orb.imag(), src, tag+10000, comm, orbinfo.nChunksImag);
    }
#endif
}

} //namespace orbital


} //namespace mrchem
