#include "SerialFunctionTree.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
#include "GenNode.h"

using namespace std;

int NFtrees=0;

/** SerialTree class constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function and allocate their nodes. 
  * NOTES:
  * Serial trees are made of projected nodes, and include gennodes and loose nodes separately.
  * All created (using class creator) Projected nodes or GenNodes are loose nodes. 
  * Loose nodes have their coeff in serial Tree, but not the node part. 
  * Projected nodes and GenNodes that are created by their creator, are detroyed by destructor ~ProjectedNode and ~GenNode. 
  * Serial tree nodes are not using the destructors, but explicitely call to deallocNodes or deallocGenNodes
  * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
*/
template<int D>
SerialFunctionTree<D>::SerialFunctionTree(FunctionTree<D> *tree, int max_nodes)
        : SerialTree<D>(tree),
          nGenNodes(0),
          maxGenNodes(max_nodes),
          lastNode(0),
          lastGenNode(0) {


    this->maxNodes = max_nodes;
    this->nNodes = 0;

    NFtrees++;
    if(MPI_rank==0 and NFtrees%10==1) println(10," N Function trees created: "<<NFtrees<<" max_nodes = " << max_nodes<<" dim = "<<D);

    //Size for GenNodes chunks. ProjectedNodes will be 8 times larger
    this->sizeGenNodeCoeff = this->tree_p->getKp1_d();//One block
    this->sizeNodeCoeff = (1<<D)*this->sizeGenNodeCoeff;//TDim  blocks
    println(10, "SizeNode Coeff (kB) " << this->sizeNodeCoeff*sizeof(double)/1024);
    println(10, "SizeGenNode Coeff (kB) " << this->sizeGenNodeCoeff*sizeof(double)/1024);

    int sizePerChunk = 1024*1024;// 1 MB small for no waisting place, but large enough so that latency and overhead work is negligible     
    if(D<3){
      //define rather from number of nodes per chunk
      this->maxNodesPerChunk = 64;
      sizePerChunk = this->maxNodesPerChunk*this->sizeNodeCoeff;
    }else{      
      this->maxNodesPerChunk = sizePerChunk/this->sizeGenNodeCoeff;
    }

    //indicate occupation of nodes
    this->nodeStackStatus = new int[this->maxNodes + 1];
    this->genNodeStackStatus = new int[this->maxGenNodes + 1];

    this->lastNode = (ProjectedNode<D>*) this->sNodes;//position of last allocated node
    this->lastGenNode = this->sGenNodes;//position of last allocated Gen node

    //initialize stacks
    for (int i = 0; i < this->maxNodes; i++) {
        this->nodeStackStatus[i] = 0;//0=unoccupied
    }
    this->nodeStackStatus[this->maxNodes] = -1;//=unavailable

    //initialize stacks
    for (int i = 0; i < this->maxGenNodes; i++) {
        this->genNodeStackStatus[i] = 0;//0=unoccupied
    }
    this->genNodeStackStatus[this->maxGenNodes] = -1;//=unavailable

    //make virtual table pointers
    ProjectedNode<D>* tmpNode = new ProjectedNode<D>();
    this->cvptr_ProjectedNode =  *(char**)(tmpNode);
    delete tmpNode;

    GenNode<D>* tmpGenNode = new GenNode<D>();
    this->cvptr_GenNode =  *(char**)(tmpGenNode);
    delete tmpGenNode;

#ifdef HAVE_OPENMP
    omp_init_lock(&Sfunc_tree_lock);
#endif
}

/** SerialTree destructor. */
template<int D>
SerialFunctionTree<D>::~SerialFunctionTree() {
    for (int i = 0; i < this->genNodeCoeffChunks.size(); i++) delete[] this->genNodeCoeffChunks[i];
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[] (char*)(this->nodeChunks[i]);
    for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    for (int i = 0; i < this->genNodeChunks.size(); i++) delete[] (char*)(this->genNodeChunks[i]);

    delete[] this->nodeStackStatus;
    delete[] this->genNodeStackStatus;

    NFtrees--;

#ifdef HAVE_OPENMP
    omp_destroy_lock(&Sfunc_tree_lock);
#endif
}

template<int D>
void SerialFunctionTree<D>::allocRoots(MWTree<D> &tree) {
    int sIx;
    double *coefs_p;
    //reserve place for nRoots
    int nRoots = tree.getRootBox().size();
    ProjectedNode<D> *root_p = this->allocNodes(nRoots, &sIx, &coefs_p);

    MWNode<D> **roots = tree.getRootBox().getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        roots[rIdx] = root_p;

        *(char**)(root_p) = this->cvptr_ProjectedNode;

        root_p->tree = &tree;
        root_p->parent = 0;
        for (int i = 0; i < root_p->getTDim(); i++) {
            root_p->children[i] = 0;
        }

        root_p->nodeIndex = tree.getRootBox().getNodeIndex(rIdx);
        root_p->hilbertPath = HilbertPath<D>();

        root_p->n_coefs = this->sizeNodeCoeff;
        root_p->coefs = coefs_p;

        root_p->serialIx = sIx;
        root_p->parentSerialIx = -1;//to indicate rootnode
        root_p->childSerialIx = -1;

        root_p->status = 0;

        root_p->clearNorms();
        root_p->setIsLeafNode();
        root_p->setIsAllocated();
        root_p->clearHasCoefs();
        root_p->setIsEndNode();
        root_p->setIsRootNode();

        tree.incrementNodeCount(root_p->getScale());

#ifdef HAVE_OPENMP
        omp_init_lock(&(root_p->node_lock));
#endif

        sIx++;
        root_p++;
        coefs_p += this->sizeNodeCoeff;
    }
}

template<int D>
void SerialFunctionTree<D>::allocChildren(MWNode<D> &parent) {
    int sIx;
    double *coefs_p;
    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    ProjectedNode<D> *child_p = this->allocNodes(nChildren, &sIx, &coefs_p);

    //position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char**)(child_p) = this->cvptr_ProjectedNode;

        child_p->tree = parent.tree;
        child_p->parent = &parent;
        for (int i = 0; i < child_p->getTDim(); i++) {
            child_p->children[i] = 0;
        }

        child_p->nodeIndex = NodeIndex<D>(parent.getNodeIndex(), cIdx);
        child_p->hilbertPath = HilbertPath<D>(parent.getHilbertPath(), cIdx);

        child_p->n_coefs = this->sizeNodeCoeff;
        child_p->coefs = coefs_p;

        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;

        child_p->status = 0;

        child_p->clearNorms();
        child_p->setIsLeafNode();
        child_p->setIsAllocated();
        child_p->clearHasCoefs();
        child_p->setIsEndNode();

        child_p->tree->incrementNodeCount(child_p->getScale());

#ifdef HAVE_OPENMP
        omp_init_lock(&child_p->node_lock);
#endif

        sIx++;
        child_p++;
        coefs_p += this->sizeNodeCoeff;
    }
}

template<int D>
void SerialFunctionTree<D>::allocGenChildren(MWNode<D> &parent) {
    int sIx;
    double *coefs_p;
    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    GenNode<D>* child_p = this->allocGenNodes(nChildren, &sIx, &coefs_p);

    //position of first child
    parent.childSerialIx = sIx;//not used fro Gennodes?
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char**)(child_p) = this->cvptr_GenNode;

	child_p->tree = parent.tree;
	child_p->parent = &parent;
	for (int i = 0; i < child_p->getTDim(); i++) {
	    child_p->children[i] = 0;
	}

	child_p->nodeIndex = NodeIndex<D>(parent.getNodeIndex(), cIdx);
	child_p->hilbertPath = HilbertPath<D>(parent.getHilbertPath(), cIdx);

	child_p->n_coefs = this->sizeGenNodeCoeff;
	child_p->coefs = coefs_p;

	child_p->serialIx = sIx;
	child_p->parentSerialIx = parent.serialIx;
	child_p->childSerialIx = -1;

	child_p->status = 0;

        child_p->clearNorms();
	child_p->setIsLeafNode();
	child_p->setIsAllocated();
	child_p->clearHasCoefs();
	child_p->setIsGenNode();

	child_p->tree->incrementGenNodeCount();

#ifdef HAVE_OPENMP
	omp_init_lock(&child_p->node_lock);
#endif

        sIx++;
        child_p++;
        coefs_p += this->sizeGenNodeCoeff;
    }
}

//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* SerialFunctionTree<D>::allocNodes(int nAlloc, int *serialIx, double **coefs_p) {
    *serialIx = this->nNodes;
    int chunkIx = *serialIx%(this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx+nAlloc > this->maxNodesPerChunk ) {
        //start on new chunk
        if (this->nNodes+nAlloc >= this->maxNodes){
            MSG_FATAL("maxNodes exceeded " << this->maxNodes);
        }

        //we want nodes allocated simultaneously to be allocated on the same pice.
        //possibly jump over the last nodes from the old chunk
        this->nNodes = this->maxNodesPerChunk*((this->nNodes+nAlloc-1)/this->maxNodesPerChunk);//start of next chunk

        int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk

        //careful: nodeChunks.size() is an unsigned int
        if (chunk+1 > this->nodeChunks.size()){
	    //need to allocate new chunk
	    this->sNodes = (ProjectedNode<D>*) new char[this->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	    this->nodeChunks.push_back(this->sNodes);
            double *sNodesCoeff = new double[this->sizeNodeCoeff*this->maxNodesPerChunk];
            this->nodeCoeffChunks.push_back(sNodesCoeff);
	    if (chunk%100==99 and D==3) println(10,endl<<" number of nodes "<<this->nNodes <<",number of Nodechunks now " << this->nodeChunks.size()<<", total size coeff  (MB) "<<(this->nNodes/1024) * this->sizeNodeCoeff/128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);
        *serialIx = this->nNodes;
        chunkIx = *serialIx%(this->maxNodesPerChunk);
    }
    assert((this->nNodes+nAlloc-1)/this->maxNodesPerChunk < this->nodeChunks.size());

    ProjectedNode<D> *newNode  = this->lastNode;
    ProjectedNode<D> *newNode_cp  = newNode;

    int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk
    *coefs_p = this->nodeCoeffChunks[chunk] + chunkIx*this->sizeNodeCoeff;
 
    for (int i = 0; i < nAlloc; i++) {
        if (this->nodeStackStatus[*serialIx+i] != 0)
	    println(0, *serialIx+i<<" NodeStackStatus: not available " << this->nodeStackStatus[*serialIx+i]);
        this->nodeStackStatus[*serialIx+i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->lastNode += nAlloc;

    return newNode;
}

template<int D>
void SerialFunctionTree<D>::deallocNodes(int serialIx) {
    if (this->nNodes < 0) {
        println(0, "minNodes exceeded " << this->nNodes);
        this->nNodes++;
    }
    this->nodeStackStatus[serialIx] = 0;//mark as available
    if (serialIx == this->nNodes-1) {//top of stack
        int topStack = this->nNodes;
        while (this->nodeStackStatus[topStack-1] == 0){
            topStack--;
            if (topStack < 1) break;
        }
        this->nNodes = topStack;//move top of stack
        //has to redefine lastNode
        int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk
        this->lastNode = this->nodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);
    }
}

//return pointer to the last active node or NULL if failed
template<int D>
GenNode<D>* SerialFunctionTree<D>::allocGenNodes(int nAlloc, int *serialIx, double **coefs_p) {
    omp_set_lock(&Sfunc_tree_lock);
    *serialIx = this->nGenNodes;
    int chunkIx = *serialIx%(this->maxNodesPerChunk);
  
    //Not necessarily wrong, but new:
    assert(nAlloc == (1<<D));

    if(chunkIx == 0  or  chunkIx+nAlloc > this->maxNodesPerChunk ){
        //start on new chunk
        if (this->nGenNodes+nAlloc >= this->maxGenNodes){
	    MSG_FATAL("maxNodes exceeded " << this->maxGenNodes);
        } 

        //we want nodes allocated simultaneously to be allocated on the same chunk.
        //possibly jump over the last nodes from the old chunk
        this->nGenNodes=this->maxNodesPerChunk*((this->nGenNodes+nAlloc-1)/this->maxNodesPerChunk);//start of next chunk

        int chunk = this->nGenNodes/this->maxNodesPerChunk;//find the right chunk

        //careful: nodeChunks.size() is an unsigned int
        if(chunk+1 > this->genNodeChunks.size()){
	    //need to allocate new chunk
	    this->sGenNodes = (GenNode<D>*) new char[this->maxNodesPerChunk*sizeof(GenNode<D>)];
	    this->genNodeChunks.push_back(this->sGenNodes);
	    double *sGenNodesCoeff = new double[this->sizeGenNodeCoeff*this->maxNodesPerChunk];
	    this->genNodeCoeffChunks.push_back(sGenNodesCoeff);
	    if(chunk%100==99 and D==3)println(10,endl<<" number of GenNodes "<<this->nGenNodes <<",number of GenNodechunks now " << this->genNodeChunks.size()<<", total size coeff  (MB) "<<(this->nGenNodes/1024) * this->sizeGenNodeCoeff/128);
        }
        this->lastGenNode = this->genNodeChunks[chunk] + this->nGenNodes%(this->maxNodesPerChunk);
        *serialIx = this->nGenNodes; 
        chunkIx = *serialIx%(this->maxNodesPerChunk);
    }
    assert((this->nGenNodes+nAlloc-1)/this->maxNodesPerChunk < this->genNodeChunks.size());

    GenNode<D>* newNode  = this->lastGenNode;
    GenNode<D>* newNode_cp  = newNode;

    int chunk = this->nGenNodes/this->maxNodesPerChunk;//find the right chunk
    *coefs_p = this->genNodeCoeffChunks[chunk] + chunkIx*this->sizeGenNodeCoeff;

    for (int i = 0; i < nAlloc; i++) {
        newNode_cp->serialIx = *serialIx+i;//Until overwritten!
        if (this->genNodeStackStatus[*serialIx+i] != 0)
	    println(0, *serialIx+i<<" NodeStackStatus: not available " << this->genNodeStackStatus[*serialIx+i]);
        this->genNodeStackStatus[*serialIx+i] = 1;
        newNode_cp++;
    }
    this->nGenNodes += nAlloc;
    this->lastGenNode += nAlloc;
    
    omp_unset_lock(&Sfunc_tree_lock);
    return newNode;
}

template<int D>
void SerialFunctionTree<D>::deallocGenNodes(int serialIx) {
    omp_set_lock(&Sfunc_tree_lock);
    if (this->nGenNodes <0) {
        println(0, "minNodes exceeded " << this->nGenNodes);
        this->nGenNodes++;
    }
    this->genNodeStackStatus[serialIx] = 0;//mark as available
    if (serialIx == this->nGenNodes-1) {//top of stack
        int topStack = this->nGenNodes;
        while (this->genNodeStackStatus[topStack-1] == 0) {
            topStack--;
            if (topStack < 1) break;
        }
        this->nGenNodes = topStack;//move top of stack
        //has to redefine lastGenNode
        int chunk = this->nGenNodes/this->maxNodesPerChunk;//find the right chunk
        this->lastGenNode = this->genNodeChunks[chunk] + this->nGenNodes%(this->maxNodesPerChunk);
    }
    omp_unset_lock(&Sfunc_tree_lock);
 }

/** Overwrite all pointers defined in the tree.
  * Necessary after sending the tree 
  * could be optimized. Should reset other counters? (GenNodes...) */
template<int D>
void SerialFunctionTree<D>::rewritePointers(int nChunks){
  //NOT_IMPLEMENTED_ABORT;
    
  int depthMax = 100;
  MWNode<D>* stack[depthMax*8];
  int slen = 0, counter = 0;

  this->nGenNodes = 0;

  //reinitialize stacks
  for (int i = 0; i < this->maxNodes; i++) {
        this->nodeStackStatus[i] = 0;
  }

  for (int i = 0; i < this->maxGenNodes; i++) {
        this->genNodeStackStatus[i] = 0;//0=unoccupied
  }
  this->genNodeStackStatus[this->maxGenNodes] = -1;//=unavailable

  this->getTree()->nNodes = 0;
  this->getTree()->nodesAtDepth.clear();
  this->getTree()->squareNorm = 0.0;

  for(int ichunk = 0 ; ichunk < nChunks; ichunk++){
    for(int inode = 0 ; inode < this->maxNodesPerChunk; inode++){
      ProjectedNode<D>* Node = (this->nodeChunks[ichunk]) + inode;
      if (Node->serialIx >= 0) {
	  //Node is part of tree, should be processed
	  this->getTree()->incrementNodeCount(Node->getScale());
	  if (Node->isEndNode()) this->getTree()->squareNorm += Node->getSquareNorm();
	  
	  //normally (intel) the virtual table does not change, but we overwrite anyway
	  *(char**)(Node) = this->cvptr_ProjectedNode;
	  
	  Node->tree = this->getTree();

	  //"adress" of coefs is the same as node, but in another array
	  Node->coefs = this->nodeCoeffChunks[ichunk]+ inode*this->sizeNodeCoeff;
	  
	  //adress of parent and children must be corrected
	  //can be on a different chunks
	  if(Node->parentSerialIx>=0){
	    int n_ichunk = Node->parentSerialIx/this->maxNodesPerChunk;
	    int n_inode = Node->parentSerialIx%this->maxNodesPerChunk;
	    Node->parent = this->nodeChunks[n_ichunk] + n_inode;
	  }else{Node->parent = 0;}
	    
	  for (int i = 0; i < Node->getNChildren(); i++) {
	    int n_ichunk = (Node->childSerialIx+i)/this->maxNodesPerChunk;
	    int n_inode = (Node->childSerialIx+i)%this->maxNodesPerChunk;
	    Node->children[i] = this->nodeChunks[n_ichunk] + n_inode;
	  }
	  this->nodeStackStatus[Node->serialIx] = 1;//occupied
#ifdef HAVE_OPENMP
	  omp_init_lock(&(Node->node_lock));
#endif
	}

    }
  }

  //update other MWTree data
  FunctionTree<D>* Tree = static_cast<FunctionTree<D>*> (this->tree_p);

  NodeBox<D> &rBox = Tree->getRootBox();
  MWNode<D> **roots = rBox.getNodes();

  for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
    roots[rIdx] = (this->nodeChunks[0]) + rIdx;//adress of roots are at start of NodeChunks[0] array
  }

  this->getTree()->resetEndNodeTable();


}


template class SerialFunctionTree<1>;
template class SerialFunctionTree<2>;
template class SerialFunctionTree<3>;
