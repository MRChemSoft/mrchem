#include "MRCPP/Printer"

#include "QMOperator.h"
#include "Orbital.h"

namespace mrchem {

void QMOperator::setApplyPrec(double prec) {
    if (this->apply_prec < 0.0) { 
        this->apply_prec = prec;
    } else if (abs(prec - this->apply_prec) > mrcpp::MachineZero) {
        MSG_ERROR("Clear operator before setup with different prec!");
    }
}

OrbitalVector QMOperator::operator()(OrbitalVector inp) {
    QMOperator &O = *this;
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital &phi_i = inp[i];
        Orbital Ophi_i = O(phi_i);
        out.push_back(Ophi_i);
    }
    return out;
}

OrbitalVector QMOperator::dagger(OrbitalVector inp) {
    QMOperator &O = *this;
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital &phi_i = inp[i];
        Orbital Ophi_i = O.dagger(phi_i);
        out.push_back(Ophi_i);
    }
    return out;
}

ComplexDouble QMOperator::operator()(Orbital bra, Orbital ket) {
    QMOperator &O = *this;
    Orbital Oket = O(ket);
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

ComplexMatrix QMOperator::operator()(OrbitalVector bra, OrbitalVector ket) {
    QMOperator &O = *this;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix result(Ni, Nj);

    result.setZero();
    for (int j = 0; j < Nj; j++) {
        Orbital &ket_j = ket[j];
        Orbital Oket_j = O(ket_j);
        for (int i = 0; i < Ni; i++) {
            Orbital &bra_i = bra[i];
            result(i, j) = orbital::dot(bra_i, Oket_j);
        }
        Oket_j.free();
    }
    return result;
}

ComplexMatrix QMOperator::dagger(OrbitalVector bra, OrbitalVector ket) {
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

ComplexDouble QMOperator::trace(OrbitalVector inp) {
    QMOperator &O = *this;
    ComplexDouble out(0.0, 0.0);
    for (int i = 0; i < inp.size(); i++) {
        Orbital &inp_i = inp[i];
        double occ_i = (double) inp_i.occ();
        out += occ_i*O(inp_i, inp_i);
    }
    return out;
}

} //namespace mrchem
