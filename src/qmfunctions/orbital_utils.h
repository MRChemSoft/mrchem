/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "mrchem.h"
#include "qmfunction_fwd.h"

namespace mrchem {
namespace orbital {

bool compare(const Orbital &phi_a, const Orbital &phi_b);
int compare_occ(const Orbital &phi_a, const Orbital &phi_b);
int compare_spin(const Orbital &phi_a, const Orbital &phi_b);

ComplexDouble dot(Orbital bra, Orbital ket);
ComplexVector dot(OrbitalVector &Bra, OrbitalVector &Ket);

void normalize(Orbital &phi);
void orthogonalize(Orbital &phi, Orbital psi);

OrbitalVector add(ComplexDouble a, OrbitalVector &Phi_a, ComplexDouble b, OrbitalVector &Phi_b, double prec = -1.0);
OrbitalVector rotate(const ComplexMatrix &U, OrbitalVector &Phi, double prec = -1.0);

OrbitalVector deep_copy(OrbitalVector &Phi);
OrbitalVector param_copy(const OrbitalVector &Phi);

OrbitalVector adjoin(OrbitalVector &Phi_a, OrbitalVector &Phi_b);
OrbitalVector disjoin(OrbitalVector &Phi, int spin);

void save_orbitals(OrbitalVector &Phi, const std::string &file, int n_orbs = -1);
OrbitalVector load_orbitals(const std::string &file, int n_orbs = -1);

void free(OrbitalVector &Phi);
void normalize(OrbitalVector &Phi);
void orthogonalize(OrbitalVector &Phi);
void orthogonalize(OrbitalVector &Phi, OrbitalVector &Psi);

ComplexMatrix calc_overlap_matrix(OrbitalVector &BraKet);
ComplexMatrix calc_overlap_matrix(OrbitalVector &Bra, OrbitalVector &Ket);
ComplexMatrix calc_lowdin_matrix(OrbitalVector &Phi);

ComplexMatrix localize(double prec, OrbitalVector &Phi);
ComplexMatrix diagonalize(double prec, OrbitalVector &Phi, ComplexMatrix &F);
ComplexMatrix orthonormalize(double prec, OrbitalVector &Phi);

int size_empty(const OrbitalVector &Phi);
int size_occupied(const OrbitalVector &Phi);
int size_singly(const OrbitalVector &Phi);
int size_doubly(const OrbitalVector &Phi);
int size_paired(const OrbitalVector &Phi);
int size_alpha(const OrbitalVector &Phi);
int size_beta(const OrbitalVector &Phi);
int get_multiplicity(const OrbitalVector &Phi);
int get_electron_number(const OrbitalVector &Phi, int spin = SPIN::Paired);

void set_spins(OrbitalVector &Phi, const IntVector &spins);
void set_errors(OrbitalVector &Phi, const DoubleVector &errors);
void set_occupancies(OrbitalVector &Phi, const IntVector &occ);

IntVector get_spins(const OrbitalVector &Phi);
IntVector get_occupancies(const OrbitalVector &Phi);
DoubleVector get_errors(const OrbitalVector &Phi);
DoubleVector get_norms(const OrbitalVector &Phi);
DoubleVector get_squared_norms(const OrbitalVector &Phi);
ComplexVector get_integrals(const OrbitalVector &Phi);

void print(const OrbitalVector &Phi);

} //namespace orbital
} //namespace mrchem
