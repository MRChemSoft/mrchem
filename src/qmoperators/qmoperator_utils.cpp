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
#include "qmoperators/QMPotential.h"
#include "qmoperators/one_electron/MomentumOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace qmoperator {
ComplexMatrix calc_kinetic_matrix_component(int d, MomentumOperator &p, OrbitalVector &bra, OrbitalVector &ket);
ComplexMatrix calc_kinetic_matrix_component(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial = false);
ComplexMatrix calc_kinetic_matrix_component_symmetrized(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial = false);
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

ComplexDouble qmoperator::calc_kinetic_trace(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &Phi, bool spinorial) {
    ComplexDouble out = {0.0, 0.0};
    int alpha_index = 0; // for spinorial operators, the kinetic operator is of the form sigma p V sigma p, with sigma being a Pauli matrix. alpha_index represents which Pauli matrix is used (x,y or z). For non-spinorial operators, it is unused.
    {   
        if (spinorial) alpha_index = 1; // if spinorial, we use sigma_x for the kinetic operator, but it is arbitrary which Pauli matrix we use, as long as it is the same on the left and right of V.
        OrbitalVector dPhi = p[0](Phi, alpha_index);
        out += V.trace(dPhi);
    }
    {   
        if (spinorial) alpha_index = 2;
        OrbitalVector dPhi = p[1](Phi, alpha_index);
        out += V.trace(dPhi);
    }
    {
        if (spinorial) alpha_index = 3;
        OrbitalVector dPhi = p[2](Phi, alpha_index);
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
ComplexMatrix qmoperator::calc_kinetic_matrix(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial) {
    ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component(0, p, V, bra, ket, spinorial); 
    ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component(1, p, V, bra, ket, spinorial);
    ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component(2, p, V, bra, ket, spinorial);
    return T_x + T_y + T_z;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_symmetrized(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial) {
    ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component_symmetrized(0, p, V, bra, ket, spinorial);
    ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component_symmetrized(1, p, V, bra, ket, spinorial);
    ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component_symmetrized(2, p, V, bra, ket, spinorial);
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
        T = mrcpp::calc_overlap_matrix(dKet);
    } else {
        OrbitalVector dBra = p[d](bra);
        OrbitalVector dKet = p[d](ket);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = mrcpp::calc_overlap_matrix(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_component_symmetrized(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);

    //in case of spinorial orbitals, the kinetic operator is sigma dot p, we therefore need to set the index of the pauli matrix to use
    int pauli_index = 0; //identity (fits non-spinorial case)
    if (spinorial) pauli_index = d + 1; //Pauli index is 1 for x, 2 for y and 3 for z, whereas the momentum operator index is 0 for x, 1 for y and 2 for z, so we need to add 1 to get the correct Pauli matrix.

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        OrbitalVector dKet = (V * p[d])(ket, pauli_index);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        T = mrcpp::calc_overlap_matrix(dKet, dKet);
    } else {
        OrbitalVector dBra = (V * p[d])(bra, pauli_index);
        OrbitalVector dKet = (V * p[d])(ket, pauli_index);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = mrcpp::calc_overlap_matrix(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_component(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);

    //in case of spinorial orbitals, the kinetic operator is sigma dot p, we therefore need to set the index of the pauli matrix to use
    int pauli_index = 0; //identity (fits non-spinorial case)
    if (spinorial) pauli_index = d + 1; //Pauli index is 1 for x, 2 for y and 3 for z, whereas the momentum operator index is 0 for x, 1 for y and 2 for z, so we need to add 1 to get the correct Pauli matrix.

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        OrbitalVector dKet = p[d](ket, pauli_index);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        T = V(dKet, dKet);
    } else {
        OrbitalVector dBra = p[d](bra, pauli_index);
        OrbitalVector dKet = p[d](ket, pauli_index);
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
