/**
 *
 *  \date July, 2016
 *  \author Peter Wind <peter.wind@uit.no> \n
 *  CTCC, University of Tromsø
 *
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <vector>

#include "parallel.h"

template<int D> class MWTree;
template<int D> class MWNode;

template<int D>
class SerialTree {
public:
    SerialTree(MWTree<D> *tree) : isShared(false), shMem(0), tree_p(tree) { }
    virtual ~SerialTree() { }

    MWTree<D>* getTree() { return this->tree_p; }

    virtual void allocRoots(MWTree<D> &tree) = 0;
    virtual void allocChildren(MWNode<D> &parent) = 0;
    virtual void allocGenChildren(MWNode<D> &parent) = 0;

    virtual void deallocNodes(int serialIx) = 0;
    virtual void deallocGenNodes(int serialIx) = 0;
    virtual void deallocGenNodeChunks() = 0;

    void S_mwTransform(double* coeff_in, double* coeff_out, bool readOnlyScaling, int stride, bool overwrite = true);
    void S_mwTransformBack(double* coeff_in, double* coeff_out, int stride);

    int nNodes;                 //number of Nodes already defined
    int maxNodesPerChunk;
    std::vector<int> nodeStackStatus;
    int sizeNodeCoeff;          //size of coeff for one node
    double **coeffStack;
    int maxNodes;               //max number of nodes that can be defined
    bool isShared;              //The coefficients are stored in shared memory
    SharedMemory *shMem;
    

protected:
    MWTree<D> *tree_p;
};

