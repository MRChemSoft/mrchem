#include "MRCPP/Printer"

#include "parallel.h"

#include "RankZeroTensorOperator.h"
#include "Orbital.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

ComplexVector RankZeroTensorOperator::getCoefVector() const {
    int nCoefs = this->coef_exp.size();
    ComplexVector out(nCoefs);
    for (int i = 0; i < nCoefs; i++) {
        out(i) = this->coef_exp[i];
    }
    return out;
}

RankZeroTensorOperator& RankZeroTensorOperator::operator=(QMOperator &O) {
    this->clear();
    this->coef_exp.push_back(1.0);
    QMOperatorVector tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

RankZeroTensorOperator& RankZeroTensorOperator::operator+=(QMOperator &O) {
    this->coef_exp.push_back(1.0);
    QMOperatorVector tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

RankZeroTensorOperator& RankZeroTensorOperator::operator-=(QMOperator &O) {
    this->coef_exp.push_back(-1.0);
    QMOperatorVector tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

RankZeroTensorOperator& RankZeroTensorOperator::operator=(const RankZeroTensorOperator &O) {
    if (this != &O) {
        this->coef_exp = O.coef_exp;
        this->oper_exp = O.oper_exp;
    }
    return *this;
}

RankZeroTensorOperator& RankZeroTensorOperator::operator+=(const RankZeroTensorOperator &O) {
    if (this != &O) {
        for (int i = 0; i < O.coef_exp.size(); i++) {
            this->coef_exp.push_back(O.coef_exp[i]);
        }
        for (int i = 0; i < O.oper_exp.size(); i++) {
            this->oper_exp.push_back(O.oper_exp[i]);
        }
    }
    return *this;
}

RankZeroTensorOperator& RankZeroTensorOperator::operator-=(const RankZeroTensorOperator &O) {
    if (this != &O) {
        for (int i = 0; i < O.coef_exp.size(); i++) {
            this->coef_exp.push_back(-O.coef_exp[i]);
        }
        for (int i = 0; i < O.oper_exp.size(); i++) {
            this->oper_exp.push_back(O.oper_exp[i]);
        }
    }
    return *this;
}

void RankZeroTensorOperator::setup(double prec) {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        for (int j = 0; j < this->oper_exp[i].size(); j++) {
            this->oper_exp[i][j]->setup(prec);
        }
    }
}

void RankZeroTensorOperator::clear() {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        for (int j = 0; j < this->oper_exp[i].size(); j++) {
            this->oper_exp[i][j]->clear();
        }
    }
}

Orbital RankZeroTensorOperator::operator()(Orbital inp) {
    if (not mpi::my_orb(inp)) return inp.paramCopy();

    OrbitalVector orb_vec;
    ComplexVector coef_vec = getCoefVector();
    for (int n = 0; n < this->oper_exp.size(); n++) {
        Orbital out_n = applyOperTerm(n, inp);
        orb_vec.push_back(out_n);
    }
    Orbital out = orbital::multiply(coef_vec, orb_vec);
    orbital::free(orb_vec);
    return out;
}

Orbital RankZeroTensorOperator::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

OrbitalVector RankZeroTensorOperator::operator()(OrbitalVector &inp) {
    RankZeroTensorOperator &O = *this;
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = O(inp[i]);
        out.push_back(out_i);
    }
    return out;
}

OrbitalVector RankZeroTensorOperator::dagger(OrbitalVector &inp) {
    NOT_IMPLEMENTED_ABORT;
}

ComplexDouble RankZeroTensorOperator::operator()(Orbital bra, Orbital ket) {
    ComplexDouble out(0.0, 0.0);
    if (mpi::my_orb(bra) and mpi::my_orb(ket)) {
        for (int n = 0; n < this->oper_exp.size(); n++) {
            Orbital Oket = applyOperTerm(n, ket);
            ComplexDouble c_n = this->coef_exp[n];
            out += c_n*orbital::dot(bra, Oket);
            Oket.free();
        }
    }
    return out;
}

ComplexDouble RankZeroTensorOperator::dagger(Orbital bra, Orbital ket) {
    NOT_IMPLEMENTED_ABORT;
}

ComplexMatrix RankZeroTensorOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    ComplexMatrix out = ComplexMatrix::Zero(bra.size(), ket.size());
    for (int n = 0; n < this->oper_exp.size(); n++) {
        OrbitalVector Oket;
        for (int j = 0; j < ket.size(); j++) {
            Orbital Oket_j = applyOperTerm(n, ket[j]);
            Oket.push_back(Oket_j);
        }
        ComplexDouble c_n = this->coef_exp[n];
        ComplexMatrix O_n = orbital::calc_overlap_matrix(bra, Oket);
        out = out + c_n*O_n;
        orbital::free(Oket);
    }
    return out;
}

ComplexMatrix RankZeroTensorOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

ComplexDouble RankZeroTensorOperator::trace(OrbitalVector &Phi) {
    RankZeroTensorOperator &O = *this;

    ComplexDouble result = 0.0;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            double eta_i = (double) Phi[i].occ();
            result += eta_i*O(Phi[i], Phi[i]);
        }
    }
#ifdef HAVE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_C_DOUBLE_COMPLEX, MPI_SUM, mpi::comm_orb);
#endif
    return result;
}

ComplexDouble RankZeroTensorOperator::trace(OrbitalVector &Phi,
                                            OrbitalVector &X,
                                            OrbitalVector &Y) {
    NEEDS_TESTING;
    RankZeroTensorOperator &O = *this;

    ComplexDouble result(0.0, 0.0);
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            if (not mpi::my_orb(X[i])) MSG_ERROR("MPI communication needed");
            if (not mpi::my_orb(Y[i])) MSG_ERROR("MPI communication needed");
            double eta_i = (double) Phi[i].occ();
            ComplexDouble result_1 = O(Phi[i], X[i]);
            ComplexDouble result_2 = O(Y[i], Phi[i]);
            result += eta_i*(result_1 + result_2);
        }
    }
#ifdef HAVE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_C_DOUBLE_COMPLEX, MPI_SUM, mpi::comm_orb);
#endif

    return result;
}

Orbital RankZeroTensorOperator::applyOperTerm(int n, Orbital inp) {
    if (n >= this->oper_exp.size()) MSG_FATAL("Invalid oper term");
    if (not mpi::my_orb(inp)) return inp.paramCopy();

    Orbital out = inp;
    for (int m = 0; m < this->oper_exp[n].size(); m++) {
        if (this->oper_exp[n][m] == 0) MSG_FATAL("Invalid oper term");
        QMOperator &O_nm = *this->oper_exp[n][m];
        Orbital tmp = O_nm.apply(out);
        if (m > 0) out.free(); // don't free input orbital
        out = tmp;
    }
    return out;
}

} //namespace mrchem
