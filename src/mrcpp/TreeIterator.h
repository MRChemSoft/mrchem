#pragma once

#include "MWTree.h"
#include "MWNode.h"
#include "constants.h"

template <int D> struct IteratorNode;

template<int D>
class TreeIterator {
public:
    TreeIterator(int dir = TopDown)
        : mode(dir),
          returnGenNodes(true),
          state(0),
          initialState(0) {
    }
    virtual ~TreeIterator() {
        if (this->initialState != 0) {
            delete this->initialState;
        }
    }

    bool next() {
        if (not this->state) {
            return false;
        }
        if (this->mode == TopDown) {
            if (this->tryNode()) {
                return true;
            }
        }
        MWNode<D> &node = *this->state->node;
        if (checkDepth(node) and checkGenerated(node)) {
            const int nChildren = 1 << D;
            for (int i = 0; i < nChildren; i++) {
                int cIdx = getChildIndex(i);
                if (this->tryChild(cIdx)) {
                    return true;
                }
            }
        }
        if (this->tryNextRoot()) {
            return true;
        }
        if (this->mode == BottomUp) {
            if (this->tryNode()) {
                return true;
            }
        }
        this->removeState();
        return next();
    }

    void setReturnGenNodes(bool i = true) { this->returnGenNodes = i; }
    void setMaxDepth(int depth) { this->maxDepth = depth; }

    MWNode<D> &getNode() { return *this->state->node; }
    friend class IteratorNode<D>;
protected:
    int root;
    int nRoots;
    int mode;
    int maxDepth;
    bool returnGenNodes;
    IteratorNode<D> *state;
    IteratorNode<D> *initialState;

    virtual int getChildIndex(int i) const = 0;

    void init(MWTree<D> *tree) {
        this->root = 0;
        this->maxDepth = -1;
        this->nRoots = tree->getRootBox().size();
        this->state = new IteratorNode<D>(&tree->getRootBox().getNode(this->root));
        // Save the first state so it can be properly deleted later
        this->initialState = this->state;
    }
    bool tryNode() {
        if (not this->state) {
            return false;
        }
        if (this->state->doneNode) {
            return false;
        }
        this->state->doneNode = true;
        return true;
    }
    bool tryChild(int i) {
        if (not this->state) {
            return false;
        }
        if (this->state->doneChild[i]) {
            return false;
        }
        this->state->doneChild[i] = true;
        if (this->state->node->isLeafNode()) {
            return false;
        }
        if (this->state->node->isLeafNode()) {
            return false;
        }
        MWNode<D> *child = &this->state->node->getMWChild(i);
        this->state = new IteratorNode<D>(child, this->state);
        return next();
    }
    bool tryNextRoot() {
        if (not this->state) {
            return false;
        }
        if (not this->state->node->isRootNode()) {
            return false;
        }
        this->root++;
        if (this->root >= this->nRoots) {
            return false;
        }
        MWNode<D> *nextRoot = &state->node->getMWTree().getRootBox().getNode(root);
        this->state = new IteratorNode<D>(nextRoot, this->state);
        return next();
    }
    void removeState() {
        if (this->state == this->initialState) {
            this->initialState = 0;
        }
        if (this->state != 0) {
            IteratorNode<D> *spare = this->state;
            this->state = spare->next;
            spare->next = 0;
            delete spare;
        }
    }
    void setDirection(int dir) {
        switch(dir) {
        case TopDown:
            this->mode = TopDown;
            break;
        case BottomUp:
            this->mode = BottomUp;
            break;
        default:
            MSG_FATAL("Invalid recursive direction!");
            break;
        }
    }
    bool checkDepth(const MWNode<D> &node) const {
        if (this->maxDepth < 0) {
            return true;
        } else if (node.getDepth() < this->maxDepth) {
            return true;
        } else {
            return false;
        }
    }
    bool checkGenerated(const MWNode<D> &node) const {
        if (node.isEndNode() and not this->returnGenNodes) {
            return false;
        } else {
            return true;
        }
    }
};

template<int D>
class IteratorNode {
public:
    MWNode<D> *node;
    IteratorNode<D> *next;
    bool doneNode;
    bool doneChild[1 << D];

    IteratorNode(MWNode<D> *nd, IteratorNode<D> *nx = 0)
            : node(nd),
              next(nx),
              doneNode(false) {
        int nChildren = 1 << D;
        for (int i = 0; i < nChildren; i++) {
            this->doneChild[i] = false;
        }
    }

    ~IteratorNode() { delete this->next; }
};

