#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "QMOperatorExp.h"
#include "Orbital.h"

template<int I, class T>
class QMTensorOperator {
public:
    QMTensorOperator() { }
    virtual ~QMTensorOperator() { }

    void setOperator(int i, T &oper) { this->oper[i] = oper; }
    void setup(double prec) { for (int i = 0; i < I; i++) this->oper[i].setup(prec); }
    void clear() { for (int i = 0; i < I; i++) this->oper[i].clear(); }

    T& operator[](int i) { return this->oper[i]; }
    const T& operator[](int i) const { return this->oper[i]; }

protected:
    T oper[I];
};

typedef QMOperatorExp RankZeroTensorOperator;

template<int I>
class RankOneTensorOperator : public QMTensorOperator<I, RankZeroTensorOperator> {
public:
    Eigen::VectorXd trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
        RankOneTensorOperator &O = *this;
        Eigen::VectorXd result = Eigen::VectorXd::Zero(I);
        for (int i = 0; i < I; i++) {
            result(i) = O[i].trace(phi, x, y);
        }
        return result;
    }
    Eigen::VectorXd trace(OrbitalVector &phi) {
        RankOneTensorOperator &O = *this;
        Eigen::VectorXd result = Eigen::VectorXd::Zero(I);
        for (int i = 0; i < I; i++) {
            result(i) = O[i].trace(phi);
        }
        return result;
    }
};

template<int I, int J>
class RankTwoTensorOperator : public QMTensorOperator<I, RankOneTensorOperator<J> > {
public:
    Eigen::MatrixXd trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
        RankTwoTensorOperator &O = *this;
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(I,J);
        for (int i = 0; i < I; i++) {
            result.row(i) = O[i].trace(phi, x, y);
        }
        return result;
    }
    Eigen::MatrixXd trace(OrbitalVector &phi) {
        RankTwoTensorOperator &O = *this;
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(I,J);
        for (int i = 0; i < I; i++) {
            result.row(i) = O[i].trace(phi);
        }
        return result;
    }
};


