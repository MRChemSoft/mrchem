/*
 * \breif Simple storage class for scale and translation indexes.
 * The usefulness of the class becomes evident when examining
 * the parallel algorithms for projection & friends.
 */

#pragma once

#include <iostream>

#include "TelePrompter.h"

template<int D> class NodeIndexComp;

template<int D>
class NodeIndex {
public:
    NodeIndex(int n = 0, const int *l = 0);
    NodeIndex(const NodeIndex<D> &idx);
    NodeIndex(const NodeIndex<D> &pIdx, int cIdx);
    virtual ~NodeIndex() { }

    inline NodeIndex<D>& operator=(const NodeIndex<D> &idx);
    inline bool operator==(const NodeIndex<D> &idx) const;
    inline bool operator!=(const NodeIndex<D> &idx) const;

    void setScale(int n) { this->N = (short int) n; }
    inline void setTranslation(const int *l);

    int getScale() const { return this->N; }
    int getTranslation(int d) const { assert(d >= 0 or d < D); return this->L[d]; }
    int *getTranslation() { return this->L; }
    const int *getTranslation() const { return this->L; }

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const NodeIndex<T> &idx);
    friend class NodeIndexComp<D>;

private:
    short int N;
    int L[D];
};

template<int D>
NodeIndex<D>::NodeIndex(int n, const int *l) {
    this->N = (short int) n;
    setTranslation(l);
}

template<int D>
NodeIndex<D>::NodeIndex(const NodeIndex<D> &idx) {
    this->N = idx.N;
    setTranslation(idx.L);
}

template<int D>
NodeIndex<D>::NodeIndex(const NodeIndex<D> &pIdx, int cIdx) {
    this->N = pIdx.N + 1;
    const int *l = pIdx.getTranslation();
    for (int d = 0; d < D; d++) {
        this->L[d] = (2 * l[d]) + ((cIdx >> d) & 1);
    }
}

template<int D>
NodeIndex<D>& NodeIndex<D>::operator=(const NodeIndex<D> &idx) {
    if (&idx == this) {
        return *this;
    }
    this->N = idx.N;
    setTranslation(idx.L);
    return *this;
}

template<int D>
void NodeIndex<D>::setTranslation(const int *l) {
    for (int d = 0; d < D; d++) {
        if (l != 0) {
            this->L[d] = l[d];
        } else {
            this->L[d] = 0;
        }
    }
}

template<int D>
bool NodeIndex<D>::operator==(const NodeIndex<D> &idx) const {
    if (this->N != idx.N) return false;
    for (int d = 0; d < D; d++) {
        if (this->L[d] != idx.L[d]) return false;
    }
    return true;
}

template<int D>
bool NodeIndex<D>::operator!=(const NodeIndex<D> &idx) const {
    if (this->N != idx.N) return true;
    for (int d = 0; d < D; d++) {
        if (this->L[d] != idx.L[d]) return true;
    }
    return false;
}

template<int D>
std::ostream& operator<<(std::ostream &o, const NodeIndex<D> &idx) {
    o << "[ " << idx.N << " | ";
    for (int d = 0; d < D - 1; d++) {
        o << idx.L[d] << ", ";
    }
    o << idx.L[D - 1] << "]";
    return o;
}


template<int D>
class NodeIndexComp {
public:
    bool operator()(const NodeIndex<D> &a, const NodeIndex<D> &b) const {
        if (a.N < b.N) {
            return true;
        }
        if (a.N > b.N) {
            return false;
        }
        for (int d = 0; d < D; d++) {
            if (a.L[d] == b.L[d]) {
                continue;
            }
            if (a.L[d] < b.L[d]) {
                return true;
            }
            return false;
        }
        return false;
    }

    bool operator()(const NodeIndex<D> *a, const NodeIndex<D> *b) const {
        if (a->N < b->N) {
            return true;
        }
        if (a->N > b->N) {
            return false;
        }
        for (int d = 0; d < D; d++) {
            if (a->L[d] == b->L[d]) {
                continue;
            }
            if (a->L[d] < b->L[d]) {
                return true;
            }
            return false;
        }
        return false;
    }
};

