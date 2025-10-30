/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <MRCPP/Parallel>

#include "RankZeroOperator.h"

#include "chemistry/Nucleus.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using QMOperator_p = std::shared_ptr<mrchem::QMOperator>;
using mrcpp::Timer;

namespace mrchem {

/** @brief return the expansion coefficients as an Eigen vector
 *
 * Converts std::vector<std::complex<double> > to Eigen::VectorXcd
 */
ComplexVector RankZeroOperator::getCoefVector() const {
    int nCoefs = this->coef_exp.size();
    ComplexVector out(nCoefs);
    for (int i = 0; i < nCoefs; i++) out(i) = this->coef_exp[i];
    return out;
}

/** @brief return the i-th term (c_i*q_i1*q_i2*...) in the expansion: sum_i c_i * prod_j q_ij */
RankZeroOperator RankZeroOperator::get(int i) {
    if (i < 0 or i >= this->size()) MSG_ABORT("Invalid operator term (i): " << i);
    RankZeroOperator out;
    out.name() = this->name();
    auto c_i = this->coef_exp[i];
    auto Q_i = this->oper_exp[i];
    out.coef_exp.push_back(c_i);
    out.oper_exp.push_back(Q_i);
    return out;
}

/** @brief return ij-the operator (q_ij) in the expansion: sum_i c_i * prod_j q_ij */
RankZeroOperator RankZeroOperator::get(int i, int j) {
    if (i < 0 or i >= this->size()) MSG_ABORT("Invalid operator term (i): " << i);
    if (j < 0 or j >= this->size(i)) MSG_ABORT("Invalid operator term (j): " << i);
    return this->oper_exp[i][j];
}

RankZeroOperator RankZeroOperator::operator()(RankZeroOperator B) {
    RankZeroOperator &A = *this;
    RankZeroOperator out;
    out.name() = A.name() + " (" + B.name() + ")";

    int a_terms = A.oper_exp.size();
    int b_terms = B.oper_exp.size();

    for (int i = 0; i < a_terms; i++) {
        ComplexDouble a_i = A.coef_exp[i];
        for (int k = 0; k < b_terms; k++) {
            ComplexDouble b_k = B.coef_exp[k];

            // we are composing (merging) only the two operators in the middle:
            // A(B) = (a_0 ... a_{J-1}) a_J(b_0) (b_1 ... b_L)
            auto &A_last = A.oper_exp[i].back();
            auto &B_first = B.oper_exp[k].front();
            QMOperatorVector AB = A_last->apply(B_first);

            QMOperatorVector tmp;
            int Na = A.oper_exp[i].size();
            int Nb = B.oper_exp[k].size();
            for (int l = 0; l < Nb - 1; l++) tmp.push_back(B.oper_exp[k][l]);
            for (auto &ab : AB) tmp.push_back(ab);
            for (int j = 1; j < Na; j++) tmp.push_back(A.oper_exp[i][j]);

            out.coef_exp.push_back(b_k * a_i);
            out.oper_exp.push_back(tmp);
        }
    }
    return out;
}

/** @brief assignment operator
 *
 * Clears the operator expansion and sets the right hand side as the only component.
 */
RankZeroOperator &RankZeroOperator::operator=(QMOperator_p O) {
    this->clear();
    this->coef_exp.push_back(1.0);
    QMOperatorVector tmp;
    tmp.push_back(O);
    this->oper_exp.push_back(tmp);
    return *this;
}

/** @brief in-place addition
 *
 * Adds a new term to the operator expansion.
 */
RankZeroOperator &RankZeroOperator::operator+=(QMOperator_p O) {
    this->coef_exp.push_back(1.0);
    QMOperatorVector tmp;
    tmp.push_back(O);
    this->oper_exp.push_back(tmp);
    return *this;
}

/** @brief in-place subtraction
 *
 * Adds a new term to the operator expansion with -1 coefficient.
 */
RankZeroOperator &RankZeroOperator::operator-=(QMOperator_p O) {
    this->coef_exp.push_back(-1.0);
    QMOperatorVector tmp;
    tmp.push_back(O);
    this->oper_exp.push_back(tmp);
    return *this;
}

/** @brief assignment operator
 *
 * This will shallow copy the operator expansion. Ownership not transferred.
 */
RankZeroOperator &RankZeroOperator::operator=(const RankZeroOperator &O) {
    if (this != &O) {
        this->name() = O.name();
        this->coef_exp = O.coef_exp;
        this->oper_exp = O.oper_exp;
    }
    return *this;
}

/** @brief in-place addition
 *
 * This will append an operator expansion. Ownership not transferred.
 */
RankZeroOperator &RankZeroOperator::operator+=(const RankZeroOperator &O) {
    if (this != &O) {
        if (this->size() == 0) {
            this->name() = O.name();
        } else {
            this->name() += " + " + O.name();
        }
        for (auto i : O.coef_exp) this->coef_exp.push_back(i);
        for (const auto &i : O.oper_exp) this->oper_exp.push_back(i);
    } else {
        MSG_ABORT("Cannot add self in place");
    }
    return *this;
}

/** @brief in-place subtraction
 *
 * This will append an operator expansion with negated coefficients.
 * Ownership not transferred.
 */
RankZeroOperator &RankZeroOperator::operator-=(const RankZeroOperator &O) {
    if (this != &O) {
        if (this->size() == 0) {
            this->name() = O.name();
        } else {
            this->name() += " - " + O.name();
        }
        for (auto i : O.coef_exp) this->coef_exp.push_back(-i);
        for (const auto &i : O.oper_exp) this->oper_exp.push_back(i);
    } else {
        MSG_ABORT("Cannot add self in place");
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
void RankZeroOperator::setup(double prec) {
    for (auto &i : this->oper_exp) {
        for (int j = 0; j < i.size(); j++) {
            i[j]->setup(prec); }
    }
}

/** @brief run clear on all operators in the expansion
 */
void RankZeroOperator::clear() {
    for (auto &i : this->oper_exp) {
        for (int j = 0; j < i.size(); j++) { i[j]->clear(); }
    }
}

ComplexDouble RankZeroOperator::operator()(const mrcpp::Coord<3> &r) const {
    const RankZeroOperator &O = *this;
    ComplexDouble out = {0.0, 0.0};
    for (int n = 0; n < O.size(); n++) {
        ComplexDouble out_n = {1.0, 0.0};
        for (auto O_nm : this->oper_exp[n]) {
            if (O_nm == nullptr) MSG_ABORT("Invalid oper term");
            out_n *= O_nm->evalf(r);
        }
        out += out_n;
    }
    return out;
}

ComplexDouble RankZeroOperator::dagger(const mrcpp::Coord<3> &r) const {
    const RankZeroOperator &O = *this;
    return std::conj(O(r));
}

/** @brief apply operator expansion to orbital
 *
 * @param inp: orbital on which to apply
 *
 * Applies each term of the operator expansion to the input orbital. First all
 * components of each term are applied consecutively, then the output of each term
 * is added upp with the corresponding coefficient.
 */
Orbital RankZeroOperator::operator()(Orbital inp) {
    if (inp.getNNodes(NUMBER::Total) == 0) return inp.paramCopy();

    RankZeroOperator &O = *this;
    std::vector<mrcpp::ComplexFunction> func_vec;
    ComplexVector coef_vec = getCoefVector();
    for (int n = 0; n < O.size(); n++) {
        Orbital out_n = O.applyOperTerm(n, inp);
        func_vec.push_back(out_n);
    }
    Orbital out = inp.paramCopy();
    mrcpp::cplxfunc::linear_combination(out, coef_vec, func_vec, -1.0);
    return out;
}

/** @brief apply the adjoint of the operator expansion to orbital
 *
 * @param inp: orbital on which to apply
 *
 * NOT IMPLEMENTED
 */
Orbital RankZeroOperator::dagger(Orbital inp) {
    if (inp.getNNodes(NUMBER::Total) == 0) return inp.paramCopy();

    RankZeroOperator &O = *this;
    std::vector<mrcpp::ComplexFunction> func_vec;
    ComplexVector coef_vec = getCoefVector();
    for (int n = 0; n < O.size(); n++) {
        Orbital out_n = O.daggerOperTerm(n, inp);
        func_vec.push_back(out_n);
    }
    Orbital out = inp.paramCopy();
    mrcpp::cplxfunc::linear_combination(out, coef_vec, func_vec, -1.0);
    return out;
}

/** @brief apply operator expansion to orbital vector
 *
 * @param inp: orbitals on which to apply
 *
 * This produces a new OrbitalVector of the same size as the input, containing
 * the corresponding output orbitals after applying the operator.
 */
OrbitalVector RankZeroOperator::operator()(OrbitalVector &inp) {
    RankZeroOperator &O = *this;
    OrbitalVector out;
    for (auto i = 0; i < inp.size(); i++) {
        Timer t1;
        Orbital out_i = O(inp[i]);
        out.push_back(out_i);
        std::stringstream o_name;
        o_name << O.name() << "|" << i << ">";
        print_utils::qmfunction(4, o_name.str(), out_i, t1);
    }
    return out;
}

/** @brief apply the adjoint of the operator expansion to orbital vector
 *
 * @param inp: orbitals on which to apply
 *
 * NOT IMPLEMENTED
 */
OrbitalVector RankZeroOperator::dagger(OrbitalVector &inp) {
    RankZeroOperator &O = *this;
    OrbitalVector out;
    for (auto i = 0; i < inp.size(); i++) {
        Timer t1;
        Orbital out_i = O.dagger(inp[i]);
        out.push_back(out_i);
        std::stringstream o_name;
        o_name << O.name() << "^dagger|" << i << ">";
        print_utils::qmfunction(4, o_name.str(), out_i, t1);
    }
    return out;
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
ComplexDouble RankZeroOperator::operator()(Orbital bra, Orbital ket) {
    RankZeroOperator &O = *this;
    Orbital Oket = O(ket);
    ComplexDouble out = orbital::dot(bra, Oket);
    return out;
}

/** @brief compute expectation value of adjoint operator
 *
 * @param bra: orbital on the bra side
 * @param ket: orbital on the ket side
 *
 * NOT IMPLEMENTED
 */
ComplexDouble RankZeroOperator::dagger(Orbital bra, Orbital ket) {
    RankZeroOperator &O = *this;
    Orbital Oket = O.dagger(ket);
    ComplexDouble out = orbital::dot(bra, Oket);
    return out;
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
ComplexMatrix RankZeroOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t1;
    RankZeroOperator &O = *this;
    OrbitalVector Oket = O(ket);
    ComplexMatrix out = orbital::calc_overlap_matrix(bra, Oket);
    std::stringstream o_name;
    o_name << "<i|" << O.name() << "|j>";
    mrcpp::print::tree(2, o_name.str(), orbital::get_n_nodes(Oket), orbital::get_size_nodes(Oket), t1.elapsed());
    return out;
}

/** @brief compute expectation matrix of adjoint operator
 *
 * @param bra: orbitals on the bra side
 * @param ket: orbitals on the ket side
 *
 * NOT IMPLEMENTED
 */
ComplexMatrix RankZeroOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t1;
    RankZeroOperator &O = *this;
    OrbitalVector Oket = O.dagger(ket);
    ComplexMatrix out = orbital::calc_overlap_matrix(bra, Oket);
    std::stringstream o_name;
    o_name << "<i|" << O.name() << "^dagger|j>";
    mrcpp::print::tree(2, o_name.str(), orbital::get_n_nodes(Oket), orbital::get_size_nodes(Oket), t1.elapsed());
    return out;
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
ComplexDouble RankZeroOperator::trace(OrbitalVector &Phi) {
    Timer t1;
    RankZeroOperator &O = *this;
    OrbitalVector OPhi = O(Phi);
    ComplexVector eta = orbital::get_occupations(Phi).cast<ComplexDouble>();
    ComplexVector phi_vec = orbital::dot(Phi, OPhi);

    std::stringstream o_name;
    o_name << "Trace " << O.name() << "(rho)";
    auto n_nodes = orbital::get_n_nodes(OPhi);
    auto n_size = orbital::get_size_nodes(OPhi);
    mrcpp::print::tree(2, o_name.str(), n_nodes, n_size, t1.elapsed());

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
ComplexDouble RankZeroOperator::trace(OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y) {
    Timer t1;
    RankZeroOperator &O = *this;

    OrbitalVector OPhi = O(Phi);
    auto y_nodes = orbital::get_n_nodes(OPhi);
    auto y_size = orbital::get_size_nodes(OPhi);
    auto y_vec = orbital::dot(Y, OPhi);
    OPhi.clear();

    OrbitalVector OX = O(X);
    auto x_nodes = orbital::get_n_nodes(OX);
    auto x_size = orbital::get_size_nodes(OX);
    auto x_vec = orbital::dot(Phi, OX);
    OX.clear();

    std::stringstream o_name;
    o_name << "Trace " << O.name() << "(rho_1)";
    mrcpp::print::tree(2, o_name.str(), std::max(x_nodes, y_nodes), std::max(x_size, y_size), t1.elapsed());

    ComplexVector eta = orbital::get_occupations(Phi).cast<ComplexDouble>();
    return eta.dot(x_vec + y_vec);
}

ComplexDouble RankZeroOperator::trace(const Nuclei &nucs) {
    Timer t1;
    RankZeroOperator &O = *this;
    ComplexVector coef_vec = getCoefVector();
    ComplexDouble out = 0.0;
    for (int n = 0; n < O.size(); n++) out += coef_vec[n] * O.traceOperTerm(n, nucs);

    std::stringstream o_name;
    o_name << "Trace " << O.name() << "(nucs)";
    mrcpp::print::tree(2, o_name.str(), 0, 0, t1.elapsed());

    return out;
}

/** @brief apply a single term of the operator expansion
 *
 * @param n: which term to apply
 * @param inp: orbital on which to apply
 *
 * This consecutively applies all components of a particular term of the operator
 * expansion to the input orbital.
 */
Orbital RankZeroOperator::applyOperTerm(int n, Orbital inp) {
    if (n >= this->oper_exp.size()) MSG_ABORT("Invalid oper term");
    if (inp.getNNodes(NUMBER::Total) == 0) return inp.paramCopy();

    Orbital out = inp;
    for (auto O_nm : this->oper_exp[n]) {
        if (O_nm == nullptr) MSG_ABORT("Invalid oper term");
        out = O_nm->apply(out);
    }
    return out;
}

Orbital RankZeroOperator::daggerOperTerm(int n, Orbital inp) {
    if (n >= this->oper_exp.size()) MSG_ABORT("Invalid oper term");
    if (inp.getNNodes(NUMBER::Total) == 0) return inp.paramCopy();

    Orbital out = inp;
    for (int i = this->oper_exp[n].size() - 1; i >= 0; i--) {
        auto O_nm = this->oper_exp[n][i];
        if (O_nm == nullptr) MSG_ABORT("Invalid oper term");
        out = O_nm->dagger(out);
    }
    return out;
}

ComplexDouble RankZeroOperator::traceOperTerm(int n, const Nuclei &nucs) {
    if (n >= this->oper_exp.size()) MSG_ABORT("Invalid oper term");

    ComplexDouble out = 0.0;
    for (const auto &nuc_k : nucs) {
        auto Z = nuc_k.getCharge();
        const auto &R = nuc_k.getCoord();
        ComplexDouble V_r = 1.0;
        for (auto O_nm : this->oper_exp[n]) {
            if (O_nm == nullptr) MSG_ABORT("Invalid oper term");
            V_r *= O_nm->evalf(R);
        }
        out += Z * V_r;
    }
    return out;
}
} // namespace mrchem
