#include "MRCPP/Printer"

#include "IdentityOperator.h"
#include "Orbital.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** Identity operator is a deep copy */
Orbital QMIdentity::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    return inp.deepCopy();
}

/** Identity operator is a deep copy */
Orbital QMIdentity::dagger(Orbital inp) {
    return apply(inp);
}

/** Overwrite default deep copy by more efficient dot product */
ComplexDouble IdentityOperator::operator()(Orbital bra, Orbital ket) {
    if (this->I.prec() < 0.0) MSG_ERROR("Uninitialized operator");
    return orbital::dot(bra, ket);
}

/** Overwrite default deep copy by more efficient dot product */
ComplexDouble IdentityOperator::dagger(Orbital bra, Orbital ket) {
    return operator()(bra, ket);
}

/** Overwrite default deep copy by calculation of overlap matrix */
ComplexMatrix IdentityOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    if (this->I.prec() < 0.0) MSG_ERROR("Uninitialized operator");

    ComplexMatrix S;
    if (&bra == &ket) {
        // In MPI it is beneficial to assume <bra| == |ket>
        S = orbital::calc_overlap_matrix(bra);
    } else {
        S = orbital::calc_overlap_matrix(bra, ket);
    }
    return S;
}

/** Overwrite default deep copy by calculation of overlap matrix */
ComplexMatrix IdentityOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    return operator()(bra, ket);
}

} //namespace mrchem
