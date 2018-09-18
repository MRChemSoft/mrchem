#pragma once

namespace mrchem {

/** @class TensorOperator
 *
 *  @brief Placeholder class for QMOperators
 *
 * This class provides a tensor notation for general QM operators. It is the base
 * for RankOneTensorOperator and RankTwoTensorOperator (can easily be extended to
 * higher ranks), and allows the components of the tensor operator to be accessed
 * in the usual way (x,y,z = integers):
 *
 * Rank 1: h[x]
 * Rank 2: h[x][y]
 * Rank 3: h[x][y][z]
 *
 * This is achieved by recursion: RankOneTensorOperator takes RankZeroTensorOperator
 * as second template, RankTwoTensorOperator takes RankOneTensorOperator as second
 * template, etc. (first argument is the dimension of the current rank). Note that
 * RankZeroTensorOperator is NOT derived from this class, but is a more general
 * operator expansion (see also @file qmoperator_fwd.h).
 *
 */

template<int I, class T>
class TensorOperator {
public:
    TensorOperator() { }
    virtual ~TensorOperator() { }

    void setup(double prec) { for (int i = 0; i < I; i++) this->oper[i].setup(prec); }
    void clear() { for (int i = 0; i < I; i++) this->oper[i].clear(); }

    T& operator[](int i) { return this->oper[i]; }
    const T& operator[](int i) const { return this->oper[i]; }

protected:
    T oper[I];
};

} //namespace mrchem
