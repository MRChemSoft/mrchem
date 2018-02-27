#pragma once

#include "MRCPP/Parallel"

namespace mrchem {

#ifdef HAVE_MPI
#include <mpi.h>
const int workOrbVecSize = 10;
#else
const int workOrbVecSize = 1;
#endif

extern int mpiOrbRank;
extern int mpiOrbSize;
extern int mpiShRank;
extern int mpiShSize;
extern int MPI_SH_group_rank;
extern int MPI_SH_group_size;

extern MPI_Comm mpiCommOrb;
extern MPI_Comm mpiCommSh;
extern MPI_Comm mpiCommSh_group;


class Orbital;
namespace orbital {

void send_orbital(Orbital &orb, int dst, int tag);
void isend_orbital(Orbital &orb, int dst, int tag, MPI_Request& request);
void recv_orbital(Orbital &orb, int src, int tag);

} //namespace orbital


class Density;
namespace density {

void send_density(Density &rho, int dst, int tag);
void isend_density(Density &rho, int dst, int tag, MPI_Request& request);
void recv_density(Density &rho, int src, int tag);

} //namespace density


} //namespace mrchem
