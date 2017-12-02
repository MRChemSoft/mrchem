#pragma once

#include <vector>
#include <complex>

#include "QMOperator.h"

class OrbitalVector;

class QMOperatorExp {
public:
    QMOperatorExp() { }
    QMOperatorExp(QMOperator &O) { *this = O; }
    QMOperatorExp(const QMOperatorExp &O) { *this = O; }
    QMOperatorExp& operator=(QMOperator &O);
    QMOperatorExp& operator=(const QMOperatorExp &O);
    QMOperatorExp& operator+=(QMOperator &O);
    QMOperatorExp& operator+=(const QMOperatorExp &O);
    virtual ~QMOperatorExp() { }

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &phi_p);
    virtual Orbital* adjoint(Orbital &phi_p);

    virtual double operator() (Orbital &phi_i, Orbital &phi_j);
    virtual double adjoint(Orbital &phi_i, Orbital &phi_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    virtual double trace(OrbitalVector &phi);
    virtual double trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y);

    friend QMOperatorExp operator*(std::complex<double> a, QMOperatorExp A);
    friend QMOperatorExp operator*(QMOperatorExp A, QMOperatorExp B);
    friend QMOperatorExp operator+(QMOperatorExp A, QMOperatorExp B);
    friend QMOperatorExp operator-(QMOperatorExp A, QMOperatorExp B);

protected:
    std::vector<std::complex<double> > coefs;
    std::vector<std::vector<QMOperator *> > oper_exp;
};

inline QMOperatorExp operator*(std::complex<double> a, QMOperatorExp A) {
    QMOperatorExp result;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        result.coefs.push_back(a*A.coefs[i]);
        result.oper_exp.push_back(A.oper_exp[i]);
    }
    return result;
}

inline QMOperatorExp operator*(QMOperatorExp A, QMOperatorExp B) {
    QMOperatorExp result;
    int a_terms = A.oper_exp.size();
    int b_terms = B.oper_exp.size();

    for (int i = 0; i < a_terms; i++) {
        std::complex<double> a_i = A.coefs[i];
        for (int k = 0; k < b_terms; k++) {
            std::complex<double> b_k = B.coefs[k];
            std::vector<QMOperator *> tmp;
            for (int l = 0; l < B.oper_exp[k].size(); l++) {
                QMOperator *B_kl = B.oper_exp[k][l];
                tmp.push_back(B_kl);
            }
            for (int j = 0; j < A.oper_exp[i].size(); j++) {
                QMOperator *A_ij = A.oper_exp[i][j];
                tmp.push_back(A_ij);
            }
            result.coefs.push_back(b_k*a_i);
            result.oper_exp.push_back(tmp);
        }
    }
    return result;
}

inline QMOperatorExp operator+(QMOperatorExp A, QMOperatorExp B) {
    QMOperatorExp result;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        result.coefs.push_back(A.coefs[i]);
        result.oper_exp.push_back(A.oper_exp[i]);
    }
    for (int i = 0; i < B.oper_exp.size(); i++) {
        result.coefs.push_back(B.coefs[i]);
        result.oper_exp.push_back(B.oper_exp[i]);
    }
    return result;
}

inline QMOperatorExp operator-(QMOperatorExp A, QMOperatorExp B) {
    QMOperatorExp result;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        result.coefs.push_back(A.coefs[i]);
        result.oper_exp.push_back(A.oper_exp[i]);
    }
    for (int i = 0; i < B.oper_exp.size(); i++) {
        result.coefs.push_back(-B.coefs[i]);
        result.oper_exp.push_back(B.oper_exp[i]);
    }
    return result;
}


