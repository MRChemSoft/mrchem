/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#pragma once

#include "qmoperators/QMOperator.h"

#include "chemistry/chemistry_fwd.h"
#include "qmfunctions/qmfunction_fwd.h"

namespace mrchem {

/** @class RankZeroOperator
 *
 *  @brief General expansion of QMOperators
 *
 * This class represents a general sum of products of fundamental QM operators,
 * which will be applied from right to left:
 *      \hat{Q} = c_1 \hat{q}_{11}\hat{q}_{12}\cdots
 *              + c_2 \hat{q}_{21}\hat{q}_{22}\cdots
 *              + \cdots
 *
 * The class has two main purposes:
 * (1)  it provides an easy, intuitive way to contruct new (more complicated) operators from
 *      simple fundamental operators by implementing the arithmetic operators (+,-,*). E.g.
 *      the x component of the angular momentum can be written explicitly from position
 *      and linear momentum as: l_x = r_y*p_z - r_z*p_y
 * (2)  it provides the public interface to all operators by implementing the operator()
 *      function, as well as expectation value and trace operation, using the (protected)
 *      apply() function of the fundamental QMOperators
 *
 * Operators can be constructed either on the fly by combining other existing operators,
 * or by defining a new derived class that contains the fundamental QMOperators (see e.g.
 * NuclearOperator which CONTAINS a NuclearPotential that is ASSIGNED to the operator itself
 * in the constructor). Note that the assignment operator is a shallow copy that does NOT
 * transfer ownership of the QMOperator pointers, so be careful when defining new operators
 * on the fly, as the underlying operators can go out of scope. This is why all implemented
 * operators contain their fundamental operators as data members.
 *
 */

class RankZeroOperator {
public:
    RankZeroOperator() {}
    RankZeroOperator(std::shared_ptr<QMOperator> O) { *this = O; }
    RankZeroOperator(const RankZeroOperator &O) { *this = O; }
    virtual ~RankZeroOperator() {}

    int size() const { return this->oper_exp.size(); }
    int size(int i) const { return this->oper_exp[i].size(); }

    std::string &name() { return this->oper_name; }
    const std::string &name() const { return this->oper_name; }

    void setup(double prec);
    void clear();

    ComplexDouble operator()(const mrcpp::Coord<3> &r) const;
    ComplexDouble dagger(const mrcpp::Coord<3> &r) const;

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
    ComplexDouble trace(const Nuclei &nucs);

    QMOperator &getRaw(int i, int j) { return *this->oper_exp[i][j]; }
    RankZeroOperator get(int i);
    RankZeroOperator get(int i, int j);

    RankZeroOperator operator()(RankZeroOperator B);

    RankZeroOperator &operator=(std::shared_ptr<QMOperator> O);
    RankZeroOperator &operator+=(std::shared_ptr<QMOperator> O);
    RankZeroOperator &operator-=(std::shared_ptr<QMOperator> O);
    RankZeroOperator &operator=(const RankZeroOperator &O);
    RankZeroOperator &operator+=(const RankZeroOperator &O);
    RankZeroOperator &operator-=(const RankZeroOperator &O);

    friend RankZeroOperator operator*(ComplexDouble a, RankZeroOperator A);
    friend RankZeroOperator operator*(RankZeroOperator A, RankZeroOperator B);
    friend RankZeroOperator operator+(RankZeroOperator A, RankZeroOperator B);
    friend RankZeroOperator operator-(RankZeroOperator A, RankZeroOperator B);

protected:
    std::string oper_name{"O"};
    std::vector<ComplexDouble> coef_exp;
    std::vector<QMOperatorVector> oper_exp;

    Orbital applyOperTerm(int n, Orbital inp);
    Orbital daggerOperTerm(int n, Orbital inp);
    ComplexDouble traceOperTerm(int n, const Nuclei &nucs);
    ComplexVector getCoefVector() const;
};

inline RankZeroOperator operator*(ComplexDouble a, RankZeroOperator A) {
    RankZeroOperator out;
    out.name() = A.name();
    for (int i = 0; i < A.oper_exp.size(); i++) {
        out.coef_exp.push_back(a * A.coef_exp[i]);
        out.oper_exp.push_back(A.oper_exp[i]);
    }
    return out;
}

inline RankZeroOperator operator*(RankZeroOperator A, RankZeroOperator B) {
    RankZeroOperator out;
    out.name() = A.name() + " * " + B.name();

    int a_terms = A.oper_exp.size();
    int b_terms = B.oper_exp.size();

    for (int i = 0; i < a_terms; i++) {
        ComplexDouble a_i = A.coef_exp[i];
        for (int k = 0; k < b_terms; k++) {
            ComplexDouble b_k = B.coef_exp[k];
            QMOperatorVector tmp;
            for (auto B_kl : B.oper_exp[k]) tmp.push_back(B_kl);
            for (auto A_ij : A.oper_exp[i]) tmp.push_back(A_ij);
            out.coef_exp.push_back(b_k * a_i);
            out.oper_exp.push_back(tmp);
        }
    }
    return out;
}

inline RankZeroOperator operator+(RankZeroOperator A, RankZeroOperator B) {
    RankZeroOperator out;
    out.name() = A.name() + " + " + B.name();
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

inline RankZeroOperator operator-(RankZeroOperator A, RankZeroOperator B) {
    RankZeroOperator out;
    out.name() = A.name() + " - " + B.name();
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

} // namespace mrchem
