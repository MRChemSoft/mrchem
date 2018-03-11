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
    for (int n = 0; n < this->oper_exp.size(); n++) {
        Orbital Oket = applyOperTerm(n, ket);
        ComplexDouble c_n = this->coef_exp[n];
        out += c_n*orbital::dot(bra, Oket);
        Oket.free();
    }
    return out;
}

ComplexDouble RankZeroTensorOperator::dagger(Orbital bra, Orbital ket) {
    NOT_IMPLEMENTED_ABORT;
}

ComplexMatrix RankZeroTensorOperator::operator()(OrbitalVector &bra,
                                                 OrbitalVector &ket) {
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix out = ComplexMatrix::Zero(Ni, Nj);

    for (int n = 0; n < this->oper_exp.size(); n++) {
        ComplexDouble c_n = this->coef_exp[n];
        for (int j = 0; j < Nj; j++) {
            Orbital Oket_j = applyOperTerm(n, ket[j]);
            for (int i = 0; i < Ni; i++) {
                out(i,j) += c_n*orbital::dot(bra[i], Oket_j);
            }
            Oket_j.free();
        }
    }
    return out;
}

ComplexMatrix RankZeroTensorOperator::dagger(OrbitalVector &bra,
                                             OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

ComplexDouble RankZeroTensorOperator::trace(OrbitalVector &phi) {
    RankZeroTensorOperator &O = *this;

    ComplexDouble result = 0.0;
    for (int i = 0; i < phi.size(); i++) {
	if (i%mpi::orb_size == mpi::orb_rank) {
	    double eta_i = (double) phi[i].occ();
	    result += eta_i*O(phi[i], phi[i]);
        }
    }
#ifdef HAVE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_C_DOUBLE_COMPLEX, MPI_SUM, mpi::comm_orb);
#endif
    return result;
}

ComplexDouble RankZeroTensorOperator::trace(OrbitalVector &phi,
                                            OrbitalVector &x,
                                            OrbitalVector &y) {
    NEEDS_TESTING;
    RankZeroTensorOperator &O = *this;

    ComplexDouble result(0.0, 0.0);
    for (int i = 0; i < phi.size(); i++) {
        if (i%mpi::orb_size == mpi::orb_rank) {
	    double eta_i = (double) phi[i].occ();
	    ComplexDouble result_1 = O(phi[i], x[i]);
	    ComplexDouble result_2 = O(y[i], phi[i]);
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

    Orbital out = inp;
    for (int m = 0; m < this->oper_exp[n].size(); m++) {
        if (this->oper_exp[n][m] == 0) MSG_FATAL("Invalid oper term");
        QMOperator &O_nm = *this->oper_exp[n][m];
        Orbital tmp = O_nm.apply(out);
        if (n > 0) out.free(); // don't free input orbital
        out = tmp;
    }
    return out;
}

} //namespace mrchem
