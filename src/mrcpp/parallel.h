#pragma once

#include "config.h"

#define EIGEN_DONT_PARALLELIZE

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_set_dynamic(n)
#define omp_set_lock(x)
#define omp_unset_lock(x)
#define omp_test_lock(x)
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

/** Share memory within a compute node
 */
class SharedMemory {
public:
    SharedMemory(int sh_size);
    ~SharedMemory();
    void allocShmem(int sh_size);
    double *sh_start_ptr; //start of shared block
    double *sh_max_ptr; //end of shared block
    double *sh_end_ptr; //end of used part
#ifdef HAVE_MPI
    MPI_Win sh_win; //MPI window object
#endif
};

#ifdef HAVE_MPI

extern MPI_Comm mpiCommOrb;
extern MPI_Comm mpiCommSh;
extern MPI_Comm mpiCommSh_group;

template<int D> class FunctionTree;

template<int D>
void Send_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm);
template<int D>
void IRcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm);
template<int D>
void ISend_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
template<int D>
void Rcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm);
void Assign_NxN(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter);
void Assign_NxN_sym(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter);

void Share_memory(int sh_size, double * &d_ptr, MPI_Win & MPI_win);

#else


#endif
