#pragma once

#include "MRCPP/Parallel"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace mrchem {

namespace omp {
extern int n_threads;
} //namespace omp

namespace mpi {
#ifdef HAVE_MPI
const int work_size = 10;
#else
const int work_size = 1;
#endif

extern int orb_rank;
extern int orb_size;
extern int share_rank;
extern int share_size;
extern int MPI_sh_group_rank;
extern int MPI_sh_group_size;

extern MPI_Comm comm_orb;
extern MPI_Comm comm_share;
extern MPI_Comm comm_sh_group;

void initialize(int argc, char **argv);
void finalize();
} //namespace mpi

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
