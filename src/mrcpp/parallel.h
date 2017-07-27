#ifndef PARALLEL_H_
#define PARALLEL_H_

#include "config.h"

template<int D> class FunctionTree;
class Orbital;

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

extern int MPI_rank;
extern int MPI_size;

void define_groups();
void MPI_Initializations();


#ifdef HAVE_MPI
#include <mpi.h>

template<int D>
void Send_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm);
template<int D>
void IRcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm);
template<int D>
void ISend_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm);
template<int D>
void Rcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm);
void Assign_NxN(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter);
void Assign_NxN_sym(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter);

//void Share_memory(MPI_Comm ncomm, MPI_Comm ncomm_sh, int sh_size, double * d_ptr);

#else


#endif
#endif /* PARALLEL_H_ */
