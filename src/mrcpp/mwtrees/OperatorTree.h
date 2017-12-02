#pragma once

#include "MWTree.h"

class OperatorTree: public MWTree<2> {
public:
    OperatorTree(const MultiResolutionAnalysis<2> &mra,
                 double np);
    virtual ~OperatorTree();

    double getNormPrecision() const { return this->normPrec; }

    void calcBandWidth(double prec = -1.0);
    void clearBandWidth();

    void setupOperNodeCache();
    void clearOperNodeCache();

    BandWidth &getBandWidth() { return *this->bandWidth; }
    const BandWidth &getBandWidth() const { return *this->bandWidth; }

    OperatorNode &getNode(int n, int l) { return *nodePtrAccess[n][l]; }
    const OperatorNode &getNode(int n, int l) const { return *nodePtrAccess[n][l]; }

    virtual void mwTransformDown(bool overwrite);
    virtual void mwTransformUp();

    friend std::ostream& operator <<(std::ostream &o, OperatorTree &tree) {
        o << std::endl << "*OperatorTree: " << tree.name << std::endl;
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

protected:
    const double normPrec;
    BandWidth *bandWidth;
    OperatorNode ***nodePtrStore;  ///< Avoids tree lookups
    OperatorNode ***nodePtrAccess; ///< Center (l=0) of node list

    void getMaxTranslations(Eigen::VectorXi &maxTransl);
};

