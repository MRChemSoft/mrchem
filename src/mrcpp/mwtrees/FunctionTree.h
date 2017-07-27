/**
*
*
*  \date Aug 14, 2009
*  \author Jonas Juselius <jonas.juselius@uit.no> \n
*  CTCC, University of Tromsø
*
*  Basic class for representing functions in a multiwavelet
* representation.
*/

#ifndef FUNCTIONTREE_H_
#define FUNCTIONTREE_H_

#include "TreeBuilder.h"
#include "MWTree.h"
#include "SerialFunctionTree.h"

class Orbital;

template<int D>
class FunctionTree: public MWTree<D> {
public:
    FunctionTree(const MultiResolutionAnalysis<D> &mra,
                 int max_nodes = MaxAllocNodes);
    virtual ~FunctionTree();

    void clear();

    double integrate() const;
    virtual double dot(const FunctionTree<D> &ket);
    virtual double evalf(const double *r);

    void getEndValues(Eigen::VectorXd &data);
    void setEndValues(Eigen::VectorXd &data);

    bool saveTree(const std::string &file);
    bool loadTree(const std::string &file);

    // In place operations
    void square();
    void power(double d);
    void normalize();
    void map(const RepresentableFunction<1> &func);

    FunctionTree<D>& operator *=(double c);
    FunctionTree<D>& operator *=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator +=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator -=(const FunctionTree<D> &tree);

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->getEndMWNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->rootBox.getNode(i)); }

    SerialFunctionTree<D>* getSerialFunctionTree() { return static_cast<SerialFunctionTree<D> *>(this->serialTree_p); }

    const FunctionNode<D> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->getEndMWNode(i)); }
    const FunctionNode<D> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->rootBox.getNode(i)); }

    template<int T>
    friend std::ostream& operator <<(std::ostream &o, FunctionTree<T> &tree);

    friend void SendRcv_Orbital(Orbital* Orb, int source, int dest, int tag);
};

template<int D>
std::ostream& operator<<(std::ostream &o, FunctionTree<D> &tree) {
    o << std::endl << "*FunctionTree: " << tree.name << std::endl;
    o << "  square norm: " << tree.squareNorm << std::endl;
    o << "  root scale: " << tree.getRootScale() << std::endl;
    o << "  order: " << tree.order << std::endl;
    o << "  nodes: " << tree.getNNodes() << std::endl;
    o << "  endNodes: " << tree.endNodeTable.size() << std::endl;
    o << "  genNodes: " << tree.getNGenNodes() << std::endl;
    o << "  nodes per scale: " << std::endl;
    for (int i = 0; i < tree.nodesAtDepth.size(); i++) {
        o << "    scale=" << i + tree.getRootScale() << "  nodes="
          << tree.nodesAtDepth[i] << std::endl;
    }
    return o;
}

#endif /* FUNCTIONTREE_H_*/
