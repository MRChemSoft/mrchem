/**
 *
 *
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#pragma once

#include "BoundingBox.h"

template<int D>
class NodeBox : public BoundingBox<D> {
public:
    NodeBox() { NOT_IMPLEMENTED_ABORT; }
    NodeBox(const NodeIndex<D> &idx, const int *nb = 0);
    NodeBox(const NodeBox<D> &box);
    NodeBox(const BoundingBox<D> &box);
    NodeBox<D> &operator=(const NodeBox<D> &box) { NOT_IMPLEMENTED_ABORT; }
    virtual ~NodeBox();

    void setNode(int idx, MWNode<D> **node);
    void clearNode(int idx) { this->nodes[idx] = 0; }

    MWNode<D> &getNode(const NodeIndex<D> &idx);
    MWNode<D> &getNode(const double *r);
    MWNode<D> &getNode(int i = 0);

    const MWNode<D> &getNode(const NodeIndex<D> &idx) const;
    const MWNode<D> &getNode(const double *r) const;
    const MWNode<D> &getNode(int i = 0) const;

    int getNOccupied() const { return this->nOccupied; }
    MWNode<D> **getNodes() { return this->nodes; }

protected:
    int nOccupied;      ///< Number of non-zero pointers in box
    MWNode<D> **nodes;  ///< Container of nodes

    void allocNodePointers();
    void deleteNodes();
};

