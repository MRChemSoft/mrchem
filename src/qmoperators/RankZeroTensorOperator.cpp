#include "MRCPP/Printer"

#include "parallel.h"

#include "RankZeroTensorOperator.h"
#include "Orbital.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief return the expansion coefficients as an Eigen vector
 *
 * Converts std::vector<std::complex<double> > to Eigen::VectorXcd
 */
ComplexVector RankZeroTensorOperator::getCoefVector() const {
    int nCoefs = this->coef_exp.size();
    ComplexVector out(nCoefs);
    for (int i = 0; i < nCoefs; i++) {
        out(i) = this->coef_exp[i];
    }
    return out;
}

/** @brief assignment operator
 *
 * Clears the operator expansion and sets the right hand side as the only component.
 */
RankZeroTensorOperator& RankZeroTensorOperator::operator=(QMOperator &O) {
    this->clear();
    this->coef_exp.push_back(1.0);
    QMOperatorVector tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

/** @brief in-place addition
 *
 * Adds a new term to the operator expansion.
 */
RankZeroTensorOperator& RankZeroTensorOperator::operator+=(QMOperator &O) {
    this->coef_exp.push_back(1.0);
    QMOperatorVector tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

/** @brief in-place subtraction
 *
 * Adds a new term to the operator expansion with -1 coefficient.
 */
RankZeroTensorOperator& RankZeroTensorOperator::operator-=(QMOperator &O) {
    this->coef_exp.push_back(-1.0);
    QMOperatorVector tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

/** @brief assignment operator
 *
 * This will shallow copy the operator expansion. Ownership not transferred.
 */
RankZeroTensorOperator& RankZeroTensorOperator::operator=(const RankZeroTensorOperator &O) {
    if (this != &O) {
        this->coef_exp = O.coef_exp;
        this->oper_exp = O.oper_exp;
    }
    return *this;
}

/** @brief in-place addition
 *
 * This will append an operator expansion. Ownership not transferred.
 */
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

/** @brief in-place subtraction
 *
 * This will append an operator expansion with negated coefficients.
 * Ownership not transferred.
 */
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

/** @brief run setup on all operators in the expansion
 *
 * @param prec: precision of application
 *
 * This prepares each of the fundamental operators in the expansion to be
 * applied with the given precision. Must be called prior to application.
 */
void RankZeroTensorOperator::setup(double prec) {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        for (int j = 0; j < this->oper_exp[i].size(); j++) {
            this->oper_exp[i][j]->setup(prec);
        }
    }
}

/** @brief run clear on all operators in the expansion
 */
void RankZeroTensorOperator::clear() {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        for (int j = 0; j < this->oper_exp[i].size(); j++) {
            this->oper_exp[i][j]->clear();
        }
    }
}

/** @brief apply operator expansion to orbital
 *
 * @param inp: orbital on which to apply
 *
 * Applies each term of the operator expansion to the input orbital. First all
 * components of each term are applied consecutively, then the output of each term
 * is added upp with the corresponding coefficient.
 */
Orbital RankZeroTensorOperator::operator()(Orbital inp) {
    if (not mpi::my_orb(inp)) return inp.paramCopy();

    OrbitalVector orb_vec;
    ComplexVector coef_vec = getCoefVector();
    for (int n = 0; n < this->oper_exp.size(); n++) {
        Orbital out_n = applyOperTerm(n, inp);
        orb_vec.push_back(out_n);
    }
    Orbital out = orbital::linear_combination(coef_vec, orb_vec);
    orbital::free(orb_vec);
    return out;
}

/** @brief apply the adjoint of the operator expansion to orbital
 *
 * @param inp: orbital on which to apply
 *
 * NOT IMPLEMENTED
 */
Orbital RankZeroTensorOperator::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief apply operator expansion to orbital vector
 *
 * @param inp: orbitals on which to apply
 *
 * This produces a new OrbitalVector of the same size as the input, containing
 * the corresponding output orbitals after applying the operator.
 */
OrbitalVector RankZeroTensorOperator::operator()(OrbitalVector &inp) {
    RankZeroTensorOperator &O = *this;
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = O(inp[i]);
        out.push_back(out_i);
    }
    return out;
}

/** @brief apply the adjoint of the operator expansion to orbital vector
 *
 * @param inp: orbitals on which to apply
 *
 * NOT IMPLEMENTED
 */
OrbitalVector RankZeroTensorOperator::dagger(OrbitalVector &inp) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief compute expectation value
 *
 * @param bra: orbital on the bra side
 * @param ket: orbital on the ket side
 *
 * This applies each term of the operator expansion separately and computes the
 * corresponding expectation value. Then the expectation values are added up with
 * the corresponding coefficient to yield the final result.
 */
ComplexDouble RankZeroTensorOperator::operator()(Orbital bra, Orbital ket) {
    RankZeroTensorOperator &O = *this;
    Orbital Oket = O(ket);
    ComplexDouble out = orbital::dot(bra, Oket);
    Oket.free();
    return out;
}

/** @brief compute expectation value of adjoint operator
 *
 * @param bra: orbital on the bra side
 * @param ket: orbital on the ket side
 *
 * NOT IMPLEMENTED
 */
ComplexDouble RankZeroTensorOperator::dagger(Orbital bra, Orbital ket) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief compute expectation matrix
 *
 * @param bra: orbitals on the bra side
 * @param ket: orbitals on the ket side
 *
 * This applies each term of the operator expansion separately to all orbitals of
 * ket vector, then computes the corresponding expectation matrix. Finally, the
 * expectation matrices are added up with the corresponding coefficient to yield
 * the final result.
 */
ComplexMatrix RankZeroTensorOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    RankZeroTensorOperator &O = *this;
    OrbitalVector Oket = O(ket);
    ComplexMatrix out = orbital::calc_overlap_matrix(bra, Oket);
    orbital::free(Oket);
    return out;
}

/** @brief compute expectation matrix of adjoint operator
 *
 * @param bra: orbitals on the bra side
 * @param ket: orbitals on the ket side
 *
 * NOT IMPLEMENTED
 */
ComplexMatrix RankZeroTensorOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief compute trace of operator expansion
 *
 * @param Phi: orbitals on which to apply
 *
 * This sums the diagonal elements of the expectation value of the operator
 * expansion applied to the density matrix:
 *      result = \sum_i n_i * <Phi_i|O|Phi_i>
 * Includes a MPI reduction operation in case of distributed orbitals.
 */
ComplexDouble RankZeroTensorOperator::trace(OrbitalVector &Phi) {
    RankZeroTensorOperator &O = *this;
    OrbitalVector OPhi = O(Phi);
    ComplexVector eta = orbital::get_occupancies(Phi).cast<ComplexDouble>();
    ComplexVector phi_vec = orbital::dot(Phi, OPhi);
    orbital::free(OPhi);
    return eta.dot(phi_vec);
}

/** @brief compute trace of operator expansion
 *
 * @param Phi: unperturbed orbitals
 * @param X: perturbed orbitals
 * @param Y: perturbed orbitals
 *
 * This sums the diagonal elements of the expectation value of the operator
 * expansion applied to the perturbed density matrix:
 *      result = \sum_i n_i * (<Phi_i|O|X_i> + <Y_i|O|Phi_i>)
 * Includes a MPI reduction operation in case of distributed orbitals.
 */
ComplexDouble RankZeroTensorOperator::trace(OrbitalVector &Phi,
                                            OrbitalVector &X,
                                            OrbitalVector &Y) {
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

/** @brief apply a single term of the operator expansion
 *
 * @param n: which term to apply
 * @param inp: orbital on which to apply
 *
 * This consecutively applies all components of a particular term of the operator
 * expansion to the input orbital.
 */
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
