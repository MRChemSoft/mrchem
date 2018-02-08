#include <iostream>
#include "parallel.h"
#include "TelePrompter.h"
#include "Timer.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
#include "SerialFunctionTree.h"
#include "SerialTree.h"
#include "InterpolatingBasis.h"
#include "MultiResolutionAnalysis.h"


using namespace std;

int mpiOrbRank = 0;
int mpiOrbSize = 1;
int mpiShRank = 0;
int mpiShSize = 1;

#ifdef HAVE_MPI
MPI_Comm mpiCommOrb;
MPI_Comm mpiCommSh;
#endif

/** sh_size in MB
 */
SharedMemory::SharedMemory(int sh_size){
    this->allocShmem(sh_size);
}

SharedMemory::~SharedMemory(){
#ifdef HAVE_MPI
    MPI_Win_free(&this->sh_win);//deallocates the memory block
#endif
}

void SharedMemory::allocShmem(int sh_size){
#ifdef HAVE_MPI
    //MPI_Aint types are used for adresses (can be larger than int).
    MPI_Aint size = (mpiShRank==0) ? 1024*1024*sh_size : 0;//rank 0 defines length of segment 
    int disp_unit = 16;//in order for the compiler to keep aligned. 
    //size is in bytes
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, mpiCommSh, &this->sh_start_ptr, &this->sh_win);
    MPI_Win_fence(0, this->sh_win);//wait until finished
    MPI_Aint qsize = 0;
    int qdisp = 0;
    MPI_Win_shared_query(this->sh_win, 0, &qsize, &qdisp, &this->sh_start_ptr);
    //    printf("me = %d, allocated: size=%zu, bytes at = %p %d\n", mpiShRank, qsize, this->sh_start_ptr, this->sh_start_ptr );
    MPI_Win_fence(0, this->sh_win);
    this->sh_max_ptr = this->sh_start_ptr + qsize/sizeof(double);
    //printf("me = %d, start = %p, max  = %p diff %d\n", mpiShRank, this->sh_start_ptr, this->sh_max_ptr, this->sh_max_ptr-this->sh_start_ptr);
    this->sh_end_ptr = this->sh_start_ptr;
    //if(mpiOrbRank==0)d_ptr[18]=2.34;
    //printf("shared: me=%d,  val=%lf\n", shrank, d_ptr[18]);
#endif
}


/** Send or receive a serial tree using MPI
 */
void MPI_Initializations(){

#ifdef HAVE_MPI
    MPI_Init(NULL, NULL);
    define_MPI_groups();
#endif
}

#ifdef HAVE_MPI
template<int D>
void Send_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  
    println(10," STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<" sending to "<<dest);
    int count = 1;
    int STreeMeta[count];
    
    if(STree->nGenNodes != 0)MSG_FATAL("sending of GenNodes not implemented");
  
    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	//set serialIx of the unused nodes to -1
	int ishift = ichunk*STree->maxNodesPerChunk;
	for (int i = 0; i < STree->maxNodesPerChunk; i++){
	    if(STree->nodeStackStatus[ishift+i] !=1 )(STree->nodeChunks[ichunk])[i].setSerialIx(-1);
	}
	MPI_Send(STree->nodeChunks[ichunk], count, MPI_BYTE, dest, tag+1+ichunk, comm);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if (not STree->isShared or not orbIsSh(dest)){
	    MPI_Send(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, dest, tag+ichunk*1000, comm);}
    }
  
    timer.stop();
    println(10," time send     " << timer); 
}
#endif

#ifdef HAVE_MPI
template<int D>
void ISend_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  
    println(10,mpiOrbRank<<" STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<" sending to "<<dest);
    int count = 1;
    int STreeMeta[count];

    timer.start();
    for (int ichunk = 0 ; ichunk < Nchunks ; ichunk++){
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	//set serialIx of the unused nodes to -1
	int ishift = ichunk*STree->maxNodesPerChunk;
	for (int i = 0; i < STree->maxNodesPerChunk; i++){
	    if (STree->nodeStackStatus[ishift+i] != 1)(STree->nodeChunks[ichunk])[i].setSerialIx(-1);
	}
	MPI_Isend(STree->nodeChunks[ichunk], count, MPI_BYTE, dest, tag+1+ichunk, comm, &request);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if (not STree->isShared or not orbIsSh(dest)){
	    MPI_Isend(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, dest, tag+ichunk*1000, comm, &request);}
    }
    timer.stop();
    println(10, " time send     " << timer); 
}
#endif

#ifdef HAVE_MPI
template<int D>
void Rcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  
    println(10, mpiOrbRank<<" STree  at "<<STree<<" receiving  "<<Nchunks<<" chunks from "<<source);
    int count = 1;
    int STreeMeta[count];

    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
	if (ichunk<STree->nodeChunks.size()) {
	    STree->sNodes = STree->nodeChunks[ichunk];
	}else{
	    double *sNodesCoeff;
	    if( STree->isShared and orbIsSh(source)) {
		if(ichunk == 0) STree->shMem->sh_end_ptr = STree->shMem->sh_start_ptr;
		sNodesCoeff = STree->shMem->sh_end_ptr;	
		STree->shMem->sh_end_ptr += (STree->sizeNodeCoeff*STree->maxNodesPerChunk);		
	    }else{
		sNodesCoeff = new double[STree->sizeNodeCoeff*STree->maxNodesPerChunk];
	    }
	    STree->nodeCoeffChunks.push_back(sNodesCoeff);
	    STree->sNodes = (ProjectedNode<D>*) new char[STree->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	    STree->nodeChunks.push_back(STree->sNodes);
	}      
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	MPI_Recv(STree->nodeChunks[ichunk], count, MPI_BYTE, source, tag+1+ichunk, comm, &status);
	println(10, mpiOrbRank<<" received  "<<count/1024<<" kB with ProjectedNodes from "<<source);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if( not STree->isShared or not orbIsSh(source)) {
	    MPI_Recv(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, source, tag+ichunk*1000, comm, &status);}
	println(10, " received  "<<count<<" coefficients from "<<source);
    }
    timer.stop();
    println(10, " time receive  " << timer);
    timer.start();
    STree->rewritePointers(Nchunks);
    timer.stop();
    println(10, " time rewrite pointers  " << timer);

}
#endif

#ifdef HAVE_MPI
template<int D>
void IRcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  

    println(10, mpiOrbRank<<" STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<" receiving from "<<source);
    int count = 1;
    int STreeMeta[count];

    MPI_Request request=MPI_REQUEST_NULL;

    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
	if (ichunk<STree->nodeChunks.size()) {
	    STree->sNodes = STree->nodeChunks[ichunk];
	}else{
	    double *sNodesCoeff = new double[STree->sizeNodeCoeff*STree->maxNodesPerChunk];
	    STree->nodeCoeffChunks.push_back(sNodesCoeff);
	    STree->sNodes = (ProjectedNode<D>*) new char[STree->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	    STree->nodeChunks.push_back(STree->sNodes);
	}      
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	MPI_Irecv(STree->nodeChunks[ichunk], count, MPI_BYTE, source, tag+1+ichunk, comm, &request);
	println(10, mpiOrbRank<<" received  "<<count/1024<<" kB with ProjectedNodes from "<<source);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if (not STree->isShared or mpiOrbRank/mpiShSize != source/mpiShSize){
	    MPI_Irecv(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, source, tag+ichunk*1000, comm, &request);}
	println(10, mpiOrbRank<<" received  "<<count<<" coefficients from "<<source);
    }


    timer.stop();
    println(10, " time receive  " << timer);
    timer.start();
    STree->rewritePointers(Nchunks);
    timer.stop();
    println(10, " time rewrite pointers  " << timer);

}
#endif


/** Assign work among MPI processes
 * For a NxN loop. Does assume symmetry
 * At each iteration data is processed, data is sent and data is received.
 * input: 
 * N, the size of one loop 
 * iter, the iteration 
 * output: 
 * which (i,j) pair to compute next
 * to whom send own data, from who receive data
 * NB: does not work for N even when N is not a multiple of mpiOrbSize
 */
void Assign_NxN_sym(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter){
    int iter = 0;
    int NBlock=N/mpiOrbSize;//number of subblocks
    int M = N%mpiOrbSize; //number of lines and column left after subblocks are taken out
    int Maxiter=-1;
    int iter_base;
    int iBlock = 0;
    int jBlock = 0;

    if(N%2==0 and N%mpiOrbSize!=0 ) MSG_FATAL("Assign_NxN_sym not implemented for this case");
    for (iter = 0; iter < *MaxIter; iter++) {
	doi[iter] = -1;
	doj[iter] = -1;
	sendto[iter] = -1;
	sendorb[iter] = -1;
	rcvorb[iter] = -1;
    }

    for (int dBlock = 0; dBlock < (N+mpiOrbSize-1)/mpiOrbSize; dBlock++) {
	//treat diagonal block first
	int imax=mpiOrbSize;//size of the diagonal block
	if(dBlock == NBlock) imax=M;
	int jmax=mpiOrbSize;
	if(dBlock == NBlock) jmax=M;
	for (int i = 0; i < imax; i++) {
	    for (int j = 0; j < jmax; j++) {

		//used in a diagonal block: 1+mpiOrbSize/2 
		//used in a block below diagonal: (1+mpiOrbSize)*0.5
		//   if even mpiOrbRank: 1+mpiOrbSize/2 if even block, and mpiOrbSize/2 if odd block, 0.5*(1+mpiOrbSize) on average, 
		//used in a block over diagonal: (1+mpiOrbSize)*0.5-1
		//number of entire blocks in a column: N/mpiOrbSize 
		//possibly used in lowest reduced block in column: M/2
		//used in all blocks in one column d (starting at 0): (1+(N/mpiOrbSize)*(1+mpiOrbSize))/2-d+M/2
		//used in all blocks up to column d, without d: d*(1+(N/mpiOrbSize)*(1+mpiOrbSize))/2-(d*(d-1))/2+d*M/2


		//(mpiOrbSize+j-i)%mpiOrbSize : for one block
		//((N+mpiOrbSize-1)/mpiOrbSize) : number of blocks in one column
		//-dBlock: The diagonal is not counted again in each block
		iter_base = (mpiOrbSize+j-i)%mpiOrbSize * ((N+mpiOrbSize-1)/mpiOrbSize)-dBlock;
		if( i == j ) iter_base = 0;
		//the last block may be smaller
		if( dBlock == NBlock )iter_base = (imax+j-i)%M * ((N+mpiOrbSize-1)/mpiOrbSize);
		//Shift for all the columns already taken into account
		iter = iter_base + dBlock*((1+((N+mpiOrbSize-1)/mpiOrbSize)*(1+mpiOrbSize))/2) - (dBlock*(dBlock-1))/2 + dBlock*M/2;

		if( (imax+j-i)%imax <= imax/2 ){
		    int i_glob = i + dBlock*mpiOrbSize;
		    int j_glob = j + dBlock*mpiOrbSize;

		    if(i == mpiOrbRank and ( imax%2!=0 or (imax+j-i)%imax < imax/2 or i < j )) {
			if(iter > Maxiter) Maxiter = iter;
			doi[iter]=i_glob;
			doj[iter]=j_glob;
			if(i!=j)rcvorb[iter]=j_glob;
			if(imax%2 == 0 and (imax+j-i)%imax == imax/2 and i < j) {
			    sendto[iter] = j;
			    sendorb[iter] = i_glob;	    
			}
		    }else if(j == mpiOrbRank and ( imax%2 != 0 or (imax+j-i)%imax < imax/2 or i < j )){
			sendto[iter] = i;
			sendorb[iter] = j_glob;
			if(imax%2 == 0 and (imax+j-i)%imax == imax/2 and i < j){
			    rcvorb[iter] = i;
			}
		    }

		    //all columns that can be generated from received data.	
		    if ((i == mpiOrbRank and (i <= j or (mpiOrbSize%2 == 0 and (imax+j-i)%imax == imax/2) ))){
			for (int icol = i; icol < N; icol+=mpiOrbSize) {
			    iBlock = icol/mpiOrbSize;
			    if ((i != j or iBlock > dBlock) and  (mpiOrbSize%2 != 0 or ((imax+j-i)%imax < imax/2) or ((iBlock+dBlock)%2 != 0 and i > j) or ((iBlock+dBlock)%2 == 0 and i < j)  )){
				if (iBlock != dBlock) {
				    int sh=1;//diagonal block is the first
				    if(iBlock>dBlock)sh=0;//diagonal block accounted for
				    if(i==j)sh=-dBlock;//only below diagonal are counted
				    if(mpiOrbSize%2 == 0 and (imax+j-i)%imax == imax/2) sh = -(iBlock+1)/2;//only half are present
				    if(mpiOrbSize%2 == 0 and (imax+j-i)%imax == imax/2 and (dBlock%2 == 0 and iBlock == 0)) sh = 1;//same iter as diagonalblock
				    if(mpiOrbSize%2 == 0 and (imax+j-i)%imax == imax/2 and (dBlock%2 != 0 and iBlock == 1)) sh = 1;//iBlock=0 same iter as diagonal block, start next iter
				    int i_glob = i + iBlock*mpiOrbSize;
				    int j_glob = j + dBlock*mpiOrbSize;
				    if ((imax+j-i)%imax <= imax/2 ) {
					doi[iter+iBlock+sh] = i_glob;
					doj[iter+iBlock+sh] = j_glob;}
				    //rcvorb[iter+iBlock+sh]=-1;//indicates "use same as before"
				    //sendorb[iter+iBlock+sh]=-1;//indicates "do not send anything"
				    //sendto[iter+iBlock+sh]=-1;//indicates "do not send anything"
				    if(iter+iBlock+sh > Maxiter) Maxiter = iter+iBlock+sh;
				}
			    }
			}
		    }else if(i == mpiOrbRank) {
			//part of blocks below blockdiagonal
			//for (jBlock = 0; jBlock < NBlock; jBlock++) {
			for (int jrow = j; jrow < N; jrow+=mpiOrbSize) {
			    jBlock = jrow/mpiOrbSize;
			    if(jBlock != dBlock  and (mpiOrbSize%2 or (imax+j-i)%imax < (imax)/2 or (jBlock%2 == 0))){
				int sh = 1;
				if(jBlock > dBlock) sh = 0;//not very elegant
				int i_glob = i + dBlock*mpiOrbSize;
				int j_glob = j + jBlock*mpiOrbSize;
				if((imax+j-i)%imax<(imax+1)/2){
				    doi[iter+jBlock+sh] = i_glob;
				    doj[iter+jBlock+sh] = j_glob;}
				//rcvorb[iter+jBlock+sh]=-1;//indicates "use same as before"
				//sendorb[iter+jBlock+sh]=-1;//indicates "do not send anything"
				//sendto[iter+jBlock+sh]=-1;//indicates "do not send anything"
				if(iter+jBlock+sh>Maxiter)Maxiter=iter+jBlock+sh;
			    }
			}
		    }
		}
	    }
	}
    }

    *MaxIter = Maxiter;
}


/** Assign work among MPI processes
 * For a NxN loop. Does not assume symmetry
 * At each iteration data is processed, data is sent and data is received.
 * input: 
 * N, the size of one loop 
 * iter, the iteration 
 * output: 
 * which (i,j) pair to compute next
 * to whom send own data, from who receive data
 */
void Assign_NxN(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter){
    int iter = 0;
    int NBlock=N/mpiOrbSize;//number of subblocks
    int M = N%mpiOrbSize; //number of lines and column left after subblocks are taken out
    int Maxiter = -1;
    int iter_base;
    int iBlock = 0;
    int jBlock = 0;

    for (iter = 0; iter < *MaxIter; iter++) {
	doi[iter] = -1;
	doj[iter] = -1;
	sendto[iter] = -1;
	sendorb[iter] = -1;
	rcvorb[iter] = -1;
    }

    for (int dBlock = 0; dBlock < (N+mpiOrbSize-1)/mpiOrbSize; dBlock++) {

	//1)  //diagonal block
	int imax=mpiOrbSize;
	if(dBlock==NBlock)imax=M;
	int jmax=mpiOrbSize;
	if(dBlock==NBlock)jmax=M;
	for (int i = 0; i < imax; i++) {
	    for (int j = 0; j < jmax; j++) {
		iter_base = (i+j)%mpiOrbSize * ((N+mpiOrbSize-1)/mpiOrbSize);
		if(dBlock==NBlock)iter_base = (i+j)%M * ((N+mpiOrbSize-1)/mpiOrbSize);
		iter = iter_base + dBlock*mpiOrbSize*((N+mpiOrbSize-1)/mpiOrbSize);
		int i_glob = i + dBlock*mpiOrbSize;
		int j_glob = j + dBlock*mpiOrbSize;
		if(i==mpiOrbRank){
		    //this processor receive orbital j_glob from MPIrank=j, and send its own orbital i_glob to MPIrank=j.
		    if(iter>Maxiter)Maxiter=iter;
		    doi[iter]=i_glob;
		    doj[iter]=j_glob;
		    if(i==j){
			sendto[iter]=-1;
			sendorb[iter]=-1;
			rcvorb[iter]=-1;
		    }else{
			sendto[iter]=j;
			sendorb[iter]=i_glob;
			rcvorb[iter]=j_glob;}// contains info both which orbital and indirectly who from
		}
		//2)	//treat all rows and columns that can be generated from received data	
		if(i<=j and i==mpiOrbRank){
		    //do the corresponding column 
		    //for (iBlock = 0; iBlock < NBlock; iBlock++) {
		    for (int icol = i; icol < N; icol+=mpiOrbSize) {
			iBlock=icol/mpiOrbSize;
			if(iBlock!=dBlock){
			    int sh=1;
			    if(iBlock>dBlock)sh=0;//not very elegant
			    int i_glob = i + iBlock*mpiOrbSize;
			    int j_glob = j + dBlock*mpiOrbSize;
			    doi[iter+iBlock+sh] = i_glob;
			    doj[iter+iBlock+sh] = j_glob;
			    rcvorb[iter+iBlock+sh]=-1;//indicates "use same as before"
			    sendorb[iter+iBlock+sh]=-1;//indicates "do not send anything"
			    sendto[iter+iBlock+sh]=-1;//indicates "do not send anything"
			    if(iter+iBlock+sh>Maxiter)Maxiter=iter+iBlock+sh;
			}
		    }
		}else if(i==mpiOrbRank){
		    //do the corresponding row
		    //for (jBlock = 0; jBlock < NBlock; jBlock++) {
		    for (int jrow = j; jrow < N; jrow+=mpiOrbSize) {
			jBlock=jrow/mpiOrbSize;
			if(jBlock!=dBlock){
			    int sh=1;
			    if(jBlock>dBlock)sh=0;//not very elegant
			    int i_glob = i + dBlock*mpiOrbSize;
			    int j_glob = j + jBlock*mpiOrbSize;
	    
			    doi[iter+jBlock+sh] = i_glob;
			    doj[iter+jBlock+sh] = j_glob;
			    rcvorb[iter+jBlock+sh]=-1;//indicates "use same as before"
			    sendorb[iter+jBlock+sh]=-1;//indicates "do not send anything"
			    sendto[iter+jBlock+sh]=-1;//indicates "do not send anything"
			    if(iter+jBlock+sh>Maxiter)Maxiter=iter+jBlock+sh;
			}
		    }
		}
	    }
	}
    }
    int NiterMax = ((N+mpiOrbSize-1)/mpiOrbSize)*((N+mpiOrbSize-1)/mpiOrbSize)*mpiOrbSize;
    if(Maxiter>NiterMax)cout<<"CHECK Assign_NxN "<<endl;
    *MaxIter = Maxiter;
}

/** Define all the different MPI groups
 */
void define_MPI_groups(){

#ifdef HAVE_MPI
    int MPI_worldsize;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_worldsize);
    println(10,"MPI world size: "<< MPI_worldsize ); 
    int MPI_worldrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_worldrank);

    //divide the world into groups
    //each group has its own group communicator definition
 
    //split world into groups that can share memory
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &mpiCommSh);

    MPI_Comm_rank(mpiCommSh, &mpiShRank);//rank in the Sh group
    MPI_Comm_size(mpiCommSh, &mpiShSize);
    println(10,"MPI SH size: "<< mpiShSize ); 

    //define a rank of the group. 
    MPI_Comm mpiCommSh_group;
    MPI_Comm_split(MPI_COMM_WORLD, mpiShRank, MPI_worldrank, &mpiCommSh_group);//mpiShRank is color (same color->in same group),  MPI_worldrank is key (orders rank within the groups)
    int MPI_SH_group_rank;
    //we define a new orbital rank, so that the orbitals within a shared memory group, have consecutive ranks.
    MPI_Comm_rank(mpiCommSh_group, &MPI_SH_group_rank);//
    mpiOrbRank = mpiShRank + MPI_SH_group_rank*MPI_worldsize;//NB: not consecutive numbers

    MPI_Comm_split(MPI_COMM_WORLD, 0, mpiOrbRank, &mpiCommOrb);//0 is color (same color->in same group), mpiOrbRank is key (orders rank in the group)
    MPI_Comm_rank(mpiCommOrb, &mpiOrbRank);
    MPI_Comm_size(mpiCommOrb, &mpiOrbSize);
    println(10,"orbital group size: "<< mpiOrbSize); 
    //        cout<<"world rank="<<MPI_worldrank<<"  Orbital rank="<<mpiOrbRank<<"  Shared mem rank="<<mpiShRank<<"  Shared memory group rank="<<MPI_SH_group_rank<<endl;

#endif
  
}
/** Tells whether an orbital is in the same shared memory group as the calling process
*   Assumes that orbital ranks within a shared memory group are consecutive
*/
bool orbIsSh(int orbRank){
    if (orbRank < mpiOrbRank-mpiShRank or  orbRank >= mpiOrbRank-mpiShRank + mpiShSize) {
	return false;
    } else { return true; }
}

#ifdef HAVE_MPI
/** Define MPI groups that share the same compute-node i.e. the same memory and array with shared memory between
 * sh_size is size in MBytes
 * returns a pointer to the shared memory*/
void Share_memory(int sh_size, double * &d_ptr, MPI_Win &win){

    //MPI_Aint types are used for adresses (can be larger than int).
    MPI_Aint size = (mpiShRank==0) ? 1024*1024*sh_size : 0;//rank 0 defines length of segment 
    //MPI_Win win = MPI_WIN_NULL;
    int disp_unit = 16;//in order for the compiler to keep aligned. 
    //size is in bytes
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, mpiCommSh, &d_ptr, &win);

    MPI_Win_fence(0, win);//wait until finished

    MPI_Aint qsize = 0;
    int qdisp = 0;
    int * qbase = NULL;
    MPI_Win_shared_query(win, 0, &qsize, &qdisp, &d_ptr);
    printf("me = %d, allocated: size=%zu, bytes at = %p %d\n", mpiShRank, qsize, d_ptr, d_ptr);
    MPI_Win_fence(0, win);
    //if(mpiOrbRank==0)d_ptr[18]=2.34;
    //printf("shared: me=%d,  val=%lf\n", shrank, d_ptr[18]);

    //MPI_Win_free(&win);//to include in destructor

} 
#endif

#ifdef HAVE_MPI
template void Rcv_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void Rcv_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void Rcv_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void Send_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int dest, int tag, MPI_Comm comm);
template void Send_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int dest, int tag, MPI_Comm comm);
template void Send_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int dest, int tag, MPI_Comm comm);
template void IRcv_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void IRcv_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void IRcv_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void ISend_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
template void ISend_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
template void ISend_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
#endif

