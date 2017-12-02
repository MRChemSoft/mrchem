/**
 *
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "parallel.h"
#include "mrcpp_declarations.h"

#include "NodeBox.h"
#include "MWNode.h"
#include "MultiResolutionAnalysis.h"
#include "SerialTree.h"

#ifdef HAVE_OPENMP
#define SET_TREE_LOCK() omp_set_lock(&this->tree_lock)
#define UNSET_TREE_LOCK() omp_unset_lock(&this->tree_lock)
#define TEST_TREE_LOCK() omp_test_lock(&this->tree_lock)
#else
#define SET_TREE_LOCK()
#define UNSET_TREE_LOCK()
#define TEST_TREE_LOCK() false
#endif

template<int D> class SerialFunctionTree;
class SerialOperatorTree;

template<int D>
class MWTree {
public:
    MWTree(const MultiResolutionAnalysis<D> &mra);
    virtual ~MWTree();

    void setZero();

    double getSquareNorm() const { return this->squareNorm; }
    void calcSquareNorm();
    void clearSquareNorm() { this->squareNorm = -1.0; }

    int getOrder() const { return this->order; }
    int getKp1() const { return this->order + 1; }
    int getKp1_d() const { return this->kp1_d; }
    int getDim() const { return D; }
    int getTDim() const { return this->tDim; }
    int getNNodes(int depth = -1) const;
    int getNEndNodes() const { return this->endNodeTable.size(); }
    int getNGenNodes();
    int getRootScale() const { return this->rootBox.getScale(); }
    int getDepth() const { return this->nodesAtDepth.size(); }

    NodeBox<D> &getRootBox() { return this->rootBox; }
    const NodeBox<D> &getRootBox() const { return this->rootBox; }
    const MultiResolutionAnalysis<D> &getMRA() const { return this->MRA; }

    void crop(double prec, double splitFac = 1.0, bool absPrec = true);
    void mwTransform(int type, bool overwrite = true);

    void setName(const std::string &n) { this->name = n; }
    const std::string &getName() const { return this->name; }

    int getRootIndex(const double *r) const { return this->rootBox.getBoxIndex(r); }
    int getRootIndex(const NodeIndex<D> &nIdx) const { return this->rootBox.getBoxIndex(nIdx); }

    MWNode<D> *findNode(const NodeIndex<D> &nIdx);
    const MWNode<D> *findNode(const NodeIndex<D> &nIdx) const;

    MWNode<D> &getNode(const NodeIndex<D> &nIdx);
    MWNode<D> &getNodeOrEndNode(const NodeIndex<D> &nIdx);
    const MWNode<D> &getNodeOrEndNode(const NodeIndex<D> &nIdx) const;

    MWNode<D> &getNode(const double *r, int depth = -1);
    MWNode<D> &getNodeOrEndNode(const double *r, int depth = -1);
    const MWNode<D> &getNodeOrEndNode(const double *r, int depth = -1) const;

    MWNode<D> &getEndMWNode(int i) { return *this->endNodeTable[i]; }
    MWNode<D> &getRootMWNode(int i) { return this->rootBox.getNode(i); }

    const MWNode<D> &getEndMWNode(int i) const { return *this->endNodeTable[i]; }
    const MWNode<D> &getRootMWNode(int i) const { return this->rootBox.getNode(i); }

    void makeNodeTable(MWNodeVector &nodeTable);
    void makeNodeTable(std::vector<MWNodeVector > &nodeTable);

    MWNodeVector* copyEndNodeTable();
    MWNodeVector* getEndNodeTable() { return &this->endNodeTable; }

    void resetEndNodeTable();
    void clearEndNodeTable() { this->endNodeTable.clear(); }

    void deleteGenerated();

    int getNThreads() const { return this->nThreads; }

    virtual bool saveTree(const std::string &file) { NOT_IMPLEMENTED_ABORT; }
    virtual bool loadTree(const std::string &file) { NOT_IMPLEMENTED_ABORT; }

    int countBranchNodes(int depth = -1);
    int countLeafNodes(int depth = -1);
    int countAllocNodes(int depth = -1);
    int countNodes(int depth = -1);
    void RecountNodes();

    SerialTree<D>* getSerialTree() { return this->serialTree_p; }

    friend class MWNode<D>;
    friend class GenNode<D>;
    friend class ProjectedNode<D>;
    friend class OperatorNode;
    friend class TreeBuilder<D>;
    friend class SerialTree<D>;
    friend class SerialFunctionTree<D>;
    friend class SerialOperatorTree;

protected:
    // Parameters that are set in construction and should never change
    const int nThreads;
    const MultiResolutionAnalysis<D> MRA;

    // Static default parameters
    const static int tDim = (1 << D);

    // Constant parameters that are derived internally
    const int order;
    const int kp1_d;

    // Parameters that are dynamic and can be set by user
    std::string name;

    SerialTree<D> *serialTree_p;

    // Tree data
    int nNodes;
    int *nGenNodes;
    double squareNorm;
    NodeBox<D> rootBox;            ///< The actual container of nodes
    MWNodeVector endNodeTable;	   ///< Final projected nodes
    std::vector<int> nodesAtDepth;  ///< Node counter

    virtual void mwTransformDown(bool overwrite);
    virtual void mwTransformUp();

    void allocNodeCounters();
    void deleteNodeCounters();

    void incrementNodeCount(int scale);
    void decrementNodeCount(int scale);
    void updateGenNodeCounts();
    void incrementGenNodeCount();
    void decrementGenNodeCount();

#ifdef HAVE_OPENMP
    omp_lock_t tree_lock;
#endif
};

