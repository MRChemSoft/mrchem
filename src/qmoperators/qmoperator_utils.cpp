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

#include "qmoperator_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "qmoperators/one_electron/MomentumOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace qmoperator {
ComplexMatrix calc_kinetic_matrix_component(int d, MomentumOperator &p, OrbitalVector &bra, OrbitalVector &ket);
ComplexMatrix calc_kinetic_matrix_component(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket);
ComplexMatrix calc_kinetic_matrix_component_symmetrized(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket);
} // namespace qmoperator

double qmoperator::calc_kinetic_trace(MomentumOperator &p, OrbitalVector &Phi) {
    DoubleVector eta = orbital::get_occupations(Phi).cast<double>();
    DoubleVector norms = DoubleVector::Zero(Phi.size());
    {
        OrbitalVector dPhi = p[0](Phi);
        norms += orbital::get_squared_norms(dPhi);
    }
    {
        OrbitalVector dPhi = p[1](Phi);
        norms += orbital::get_squared_norms(dPhi);
    }
    {
        OrbitalVector dPhi = p[2](Phi);
        norms += orbital::get_squared_norms(dPhi);
    }
    return 0.5 * eta.dot(norms);
}

ComplexDouble qmoperator::calc_kinetic_trace(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &Phi) {
    ComplexDouble out = {0.0, 0.0};
    {
        OrbitalVector dPhi = p[0](Phi);
        out += V.trace(dPhi);
    }
    {
        OrbitalVector dPhi = p[1](Phi);
        out += V.trace(dPhi);
    }
    {
        OrbitalVector dPhi = p[2](Phi);
        out += V.trace(dPhi);
    }
    return 0.5 * out;
}

/** @brief Expectation value matrix: T_ij = <i|T|j> = <i|p p|j>
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * Instead of applying the full kinetic operator on the ket's, the momentum
 * operator is applied both to the left and right, thus taking advantage
 * of symmetry and getting away with only first-derivative operators.
 */
ComplexMatrix qmoperator::calc_kinetic_matrix(MomentumOperator &p, OrbitalVector &bra, OrbitalVector &ket) {
    ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component(0, p, bra, ket);
    ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component(1, p, bra, ket);
    ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component(2, p, bra, ket);
    return T_x + T_y + T_z;
}

/** @brief Expectation value matrix ZORA: T_ij = <i|T_zora|j> = <i|p V_zora p|j>
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * Instead of applying the full kinetic operator on the ket's, the momentum
 * operator is applied both to the left and right, thus taking advantage
 * of symmetry and getting away with only first-derivative operators.
 */
ComplexMatrix qmoperator::calc_kinetic_matrix(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket) {
    ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component(0, p, V, bra, ket);
    ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component(1, p, V, bra, ket);
    ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component(2, p, V, bra, ket);
    return T_x + T_y + T_z;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_symmetrized(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket) {
    ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component_symmetrized(0, p, V, bra, ket);
    ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component_symmetrized(1, p, V, bra, ket);
    ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component_symmetrized(2, p, V, bra, ket);
    return T_x + T_y + T_z;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_component(int d, MomentumOperator &p, OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        OrbitalVector dKet = p[d](ket);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        T = orbital::calc_overlap_matrix(dKet);
    } else {
        OrbitalVector dBra = p[d](bra);
        OrbitalVector dKet = p[d](ket);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = orbital::calc_overlap_matrix(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_component_symmetrized(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        OrbitalVector dKet = (V * p[d])(ket);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        T = orbital::calc_overlap_matrix(dKet, dKet);
    } else {
        OrbitalVector dBra = (V * p[d])(bra);
        OrbitalVector dKet = (V * p[d])(ket);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = orbital::calc_overlap_matrix(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_component(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        OrbitalVector dKet = p[d](ket);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        T = V(dKet, dKet);
    } else {
        OrbitalVector dBra = p[d](bra);
        OrbitalVector dKet = p[d](ket);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = V(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

} // namespace mrchem
