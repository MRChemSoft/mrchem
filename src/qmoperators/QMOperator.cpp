#include "MRCPP/Printer"

#include "QMOperator.h"
#include "Orbital.h"

namespace mrchem {

void QMOperator::setApplyPrec(double prec) {
    if (this->apply_prec < 0.0) { 
        this->apply_prec = prec;
    } else if (not isSetup(prec)) {
        MSG_ERROR("Clear operator before setup with different prec!");
    }
}

bool QMOperator::isSetup(double prec) {
    double dPrec = std::abs(this->apply_prec - prec);
    double thrs = mrcpp::MachineZero;
    return (dPrec < thrs) ? true : false;
}

ComplexDouble QMOperator::apply(Orbital bra, Orbital ket) {
    QMOperator &O = *this;
    Orbital Oket = O.apply(ket);
    ComplexDouble result = orbital::dot(bra, Oket);
    Oket.free();
    return result;
}

ComplexDouble QMOperator::dagger(Orbital bra, Orbital ket) {
    QMOperator &O = *this;
    Orbital Oket = O.dagger(ket);
    ComplexDouble result = orbital::dot(bra, Oket);
    Oket.free();
    return result;
}

ComplexMatrix QMOperator::apply(OrbitalVector &bra, OrbitalVector &ket) {
    QMOperator &O = *this;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix result(Ni, Nj);
    result.setZero();
    for (int j = 0; j < Nj; j++) {
        Orbital &ket_j = ket[j];
        Orbital Oket_j = O.apply(ket_j);
        for (int i = 0; i < Ni; i++) {
            Orbital &bra_i = bra[i];
            result(i, j) = orbital::dot(bra_i, Oket_j);
        }
        Oket_j.free();
    }
    return result;
}

ComplexMatrix QMOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NEEDS_TESTING;
    QMOperator &O = *this;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix result(Ni, Nj);

    result.setZero();
    for (int j = 0; j < Nj; j++) {
        Orbital &ket_j = ket[j];
        Orbital Oket_j = O.dagger(ket_j);
        for (int i = 0; i < Ni; i++) {
            Orbital &bra_i = bra[i];
            result(i, j) = orbital::dot(bra_i, Oket_j);
        }
        Oket_j.free();
    }
    return result;
}

} //namespace mrchem
