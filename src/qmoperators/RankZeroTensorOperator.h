#pragma once

#include "QMOperator.h"

namespace mrchem {

// Convenience typedef
typedef std::vector<QMOperator *> QMOperatorVector;

class RankZeroTensorOperator {
public:
    RankZeroTensorOperator() { }
    RankZeroTensorOperator(QMOperator &O) { *this = O; }
    RankZeroTensorOperator(const RankZeroTensorOperator &O) { *this = O; }
    virtual ~RankZeroTensorOperator() { }

    void setup(double prec);
    void clear();

    Orbital operator()(Orbital inp);
    Orbital dagger(Orbital inp);

    OrbitalVector operator()(OrbitalVector &inp);
    OrbitalVector dagger(OrbitalVector &inp);

    ComplexDouble operator()(Orbital bra, Orbital ket);
    ComplexDouble dagger(Orbital bra, Orbital ket);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);

    ComplexDouble trace(OrbitalVector &Phi);
    ComplexDouble trace(OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y);

    RankZeroTensorOperator& operator=(QMOperator &O);
    RankZeroTensorOperator& operator+=(QMOperator &O);
    RankZeroTensorOperator& operator-=(QMOperator &O);
    RankZeroTensorOperator& operator=(const RankZeroTensorOperator &O);
    RankZeroTensorOperator& operator+=(const RankZeroTensorOperator &O);
    RankZeroTensorOperator& operator-=(const RankZeroTensorOperator &O);

    friend RankZeroTensorOperator operator*(ComplexDouble a, RankZeroTensorOperator A);
    friend RankZeroTensorOperator operator*(RankZeroTensorOperator A, RankZeroTensorOperator B);
    friend RankZeroTensorOperator operator+(RankZeroTensorOperator A, RankZeroTensorOperator B);
    friend RankZeroTensorOperator operator-(RankZeroTensorOperator A, RankZeroTensorOperator B);

protected:
    std::vector<ComplexDouble> coef_exp;
    std::vector<QMOperatorVector> oper_exp;

    Orbital applyOperTerm(int n, Orbital inp);
    ComplexVector getCoefVector() const;
};


inline RankZeroTensorOperator operator*(ComplexDouble a, RankZeroTensorOperator A) {
    RankZeroTensorOperator out;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        out.coef_exp.push_back(a*A.coef_exp[i]);
        out.oper_exp.push_back(A.oper_exp[i]);
    }
    return out;
}

inline RankZeroTensorOperator operator*(RankZeroTensorOperator A,
                                        RankZeroTensorOperator B) {
    RankZeroTensorOperator out;
    int a_terms = A.oper_exp.size();
    int b_terms = B.oper_exp.size();

    for (int i = 0; i < a_terms; i++) {
        ComplexDouble a_i = A.coef_exp[i];
        for (int k = 0; k < b_terms; k++) {
            ComplexDouble b_k = B.coef_exp[k];
            QMOperatorVector tmp;
            for (int l = 0; l < B.oper_exp[k].size(); l++) {
                QMOperator *B_kl = B.oper_exp[k][l];
                tmp.push_back(B_kl);
            }
            for (int j = 0; j < A.oper_exp[i].size(); j++) {
                QMOperator *A_ij = A.oper_exp[i][j];
                tmp.push_back(A_ij);
            }
            out.coef_exp.push_back(b_k*a_i);
            out.oper_exp.push_back(tmp);
        }
    }
    return out;
}

inline RankZeroTensorOperator operator+(RankZeroTensorOperator A,
                                        RankZeroTensorOperator B) {
    RankZeroTensorOperator out;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        out.coef_exp.push_back(A.coef_exp[i]);
        out.oper_exp.push_back(A.oper_exp[i]);
    }
    for (int i = 0; i < B.oper_exp.size(); i++) {
        out.coef_exp.push_back(B.coef_exp[i]);
        out.oper_exp.push_back(B.oper_exp[i]);
    }
    return out;
}

inline RankZeroTensorOperator operator-(RankZeroTensorOperator A,
                                        RankZeroTensorOperator B) {
    RankZeroTensorOperator out;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        out.coef_exp.push_back(A.coef_exp[i]);
        out.oper_exp.push_back(A.oper_exp[i]);
    }
    for (int i = 0; i < B.oper_exp.size(); i++) {
        out.coef_exp.push_back(-B.coef_exp[i]);
        out.oper_exp.push_back(B.oper_exp[i]);
    }
    return out;
}

} //namespace mrchem
