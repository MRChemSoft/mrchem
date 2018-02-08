#pragma once

#include "config.h"

#define EIGEN_DONT_PARALLELIZE

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


void MPI_Initializations();
void define_MPI_groups();
bool orbIsSh(int orbRank);

#ifdef HAVE_MPI

extern MPI_Comm mpiCommOrb;
extern MPI_Comm mpiCommSh;
extern MPI_Comm mpiCommSh_group;

#else


#endif
