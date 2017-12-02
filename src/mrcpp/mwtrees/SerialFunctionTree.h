/**
*
*
*  \date Jul, 2016
*  \author Peter Wind <peter.wind@uit.no> \n
*  CTCC, University of Tromsø
*
*/

#pragma once

#include <vector>

#include "SerialTree.h"
#include "parallel.h"

template<int D> class FunctionTree;
template<int D> class ProjectedNode;
template<int D> class GenNode;

template<int D>
class SerialFunctionTree : public SerialTree<D> {
public:
    SerialFunctionTree(FunctionTree<D> *tree);
    virtual ~SerialFunctionTree();

    virtual void allocRoots(MWTree<D> &tree);
    virtual void allocChildren(MWNode<D> &parent);
    virtual void allocGenChildren(MWNode<D> &parent);

    virtual void deallocNodes(int serialIx);
    virtual void deallocGenNodes(int serialIx);
    virtual void deallocGenNodeChunks();

    std::vector<ProjectedNode<D>*> nodeChunks;
    std::vector<double*> nodeCoeffChunks;

    ProjectedNode<D> *sNodes;   //serial ProjectedNodes
    GenNode<D> *sGenNodes;      //serial GenNodes

    std::vector<GenNode<D>*> genNodeChunks;
    std::vector<double*> genNodeCoeffChunks;

    int nGenNodes;              //number of GenNodes already defined

    double **genCoeffStack;

    //    int *genNodeStackStatus;
    std::vector<int> genNodeStackStatus;

    void rewritePointers(int nChunks);

protected:
    int maxGenNodes;            //max number of Gen nodes that can be defined
    int sizeGenNodeCoeff;       //size of coeff for one Gen node

    char *cvptr_ProjectedNode;  //virtual table pointer for ProjectedNode
    char *cvptr_GenNode;        //virtual table pointer for GenNode

    ProjectedNode<D>* lastNode; //pointer to the last active node
    GenNode<D>* lastGenNode;    //pointer to the last active Gen node

    ProjectedNode<D>* allocNodes(int nAlloc, int* serialIx, double **coefs_p);
    GenNode<D>* allocGenNodes(int nAlloc, int* serialIx, double **coefs_p);

private:
#ifdef HAVE_OPENMP
    omp_lock_t Sfunc_tree_lock;
#endif
};

