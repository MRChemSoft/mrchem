#include "SerialOperatorTree.h"
#include "OperatorTree.h"
#include "OperatorNode.h"

using namespace std;

int NOtrees=0;

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
SerialOperatorTree::SerialOperatorTree(OperatorTree *tree, int max_nodes)
        : SerialTree<2>(tree),
          lastNode(0) {

    this->maxNodes = max_nodes;
    this->nNodes = 0;
    NOtrees++;

    this->sizeNodeCoeff = 4*this->tree_p->getKp1_d();

    this->maxNodesPerChunk = 64;
    int sizePerChunk = this->maxNodesPerChunk*this->sizeNodeCoeff;
   
    if(MPI_rank==0 and NOtrees%100==1)println(10, " max_nodes = " << max_nodes << ", nodes per chunk = " << this->maxNodesPerChunk<<" sizePerChunk "<<sizePerChunk<<" N Op trees: "<<NOtrees);

    //indicate occupation of nodes
    this->nodeStackStatus = new int[this->maxNodes + 1];

    this->lastNode = (OperatorNode*) this->sNodes;//position of last allocated node

    //initialize stacks
    for (int i = 0; i < maxNodes; i++) {
        this->nodeStackStatus[i] = 0;//0=unoccupied
    }
    this->nodeStackStatus[maxNodes] = -1;//=unavailable

    //make virtual table pointers
    OperatorNode* tmpNode = new OperatorNode();
    this->cvptr_OperatorNode =  *(char**)(tmpNode);
    delete tmpNode;

#ifdef HAVE_OPENMP
    omp_init_lock(&Soper_tree_lock);
#endif
}

/** SerialTree destructor. */
SerialOperatorTree::~SerialOperatorTree() {
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[] (char*)(this->nodeChunks[i]);
    for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];

    delete[] this->nodeStackStatus;
    NOtrees--;

#ifdef HAVE_OPENMP
    omp_destroy_lock(&Soper_tree_lock);
#endif
}

void SerialOperatorTree::allocRoots(MWTree<2> &tree) {
    int sIx;
    double *coefs_p;
    //reserve place for nRoots
    int nRoots = tree.getRootBox().size();
    OperatorNode *root_p = this->allocNodes(nRoots, &sIx, &coefs_p);

    MWNode<2> **roots = tree.getRootBox().getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        roots[rIdx] = root_p;

        *(char**)(root_p) = this->cvptr_OperatorNode;

        root_p->tree = &tree;
        root_p->parent = 0;
        for (int i = 0; i < root_p->getTDim(); i++) {
            root_p->children[i] = 0;
        }

        root_p->nodeIndex = tree.getRootBox().getNodeIndex(rIdx);
        root_p->hilbertPath = HilbertPath<2>();

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

void SerialOperatorTree::allocChildren(MWNode<2> &parent) {
    int sIx;
    double *coefs_p;
    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    OperatorNode *child_p = this->allocNodes(nChildren, &sIx, &coefs_p);

    //position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char**)(child_p) = this->cvptr_OperatorNode;

        child_p->tree = parent.tree;
        child_p->parent = &parent;
        for (int i = 0; i < child_p->getTDim(); i++) {
            child_p->children[i] = 0;
        }

        child_p->nodeIndex = NodeIndex<2>(parent.getNodeIndex(), cIdx);
        child_p->hilbertPath = HilbertPath<2>(parent.getHilbertPath(), cIdx);

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

void SerialOperatorTree::allocGenChildren(MWNode<2> &parent) {
    NOT_REACHED_ABORT;
}

//return pointer to the last active node or NULL if failed
OperatorNode* SerialOperatorTree::allocNodes(int nAlloc, int *serialIx, double **coefs_p) {
    *serialIx = this->nNodes;
    int chunkIx = (*serialIx)%(this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx+nAlloc > this->maxNodesPerChunk ) {
        //start on new chunk
        if (this->nNodes+nAlloc >= this->maxNodes){
            MSG_FATAL("maxNodes exceeded " << this->maxNodes);
        }

        //we want nodes allocated simultaneously to be allocated on the same piece.
        //possibly jump over the last nodes from the old chunk
        this->nNodes = this->maxNodesPerChunk*((this->nNodes+nAlloc-1)/this->maxNodesPerChunk);//start of next chunk

        int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk

        //careful: nodeChunks.size() is an unsigned int
        if (chunk+1 > this->nodeChunks.size()){
	    //need to allocate new chunk
	    this->sNodes = (OperatorNode*) new char[this->maxNodesPerChunk*sizeof(OperatorNode)];
	    this->nodeChunks.push_back(this->sNodes);
            double *sNodesCoeff = new double[this->sizeNodeCoeff*this->maxNodesPerChunk];
            this->nodeCoeffChunks.push_back(sNodesCoeff);
        }
        this->lastNode = this->nodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);
        *serialIx = this->nNodes;
        chunkIx = *serialIx%(this->maxNodesPerChunk);
    }
    assert((this->nNodes+nAlloc-1)/this->maxNodesPerChunk < this->nodeChunks.size());

    OperatorNode *newNode  = this->lastNode;
    OperatorNode *newNode_cp  = newNode;

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

void SerialOperatorTree::deallocNodes(int serialIx) {
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

void SerialOperatorTree::deallocGenNodes(int serialIx) {
    NOT_REACHED_ABORT;
}
