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

ComplexDouble dot(Orbital bra, Orbital ket);
ComplexVector dot(OrbitalVector &bra, OrbitalVector &ket);

bool compare(const Orbital &orb_a, const Orbital &orb_b);
int compare_occ(const Orbital &orb_a, const Orbital &orb_b);
int compare_spin(const Orbital &orb_a, const Orbital &orb_b);

OrbitalVector add(ComplexDouble a, OrbitalVector &inp_a, ComplexDouble b, OrbitalVector &inp_b, double prec = -1.0);
OrbitalVector rotate(const ComplexMatrix &U, OrbitalVector &inp, double prec = -1.0);

OrbitalVector deep_copy(OrbitalVector &inp);
OrbitalVector param_copy(const OrbitalVector &inp);

OrbitalVector adjoin(OrbitalVector &inp_a, OrbitalVector &inp_b);
OrbitalVector disjoin(OrbitalVector &inp, int spin);

void save_orbitals(OrbitalVector &Phi, const std::string &file, int n_orbs = -1);
OrbitalVector load_orbitals(const std::string &file, int n_orbs = -1);

void free(OrbitalVector &vec);
void normalize(OrbitalVector &vec);
void orthogonalize(OrbitalVector &vec);
void orthogonalize(OrbitalVector &vec, OrbitalVector &inp);

ComplexMatrix calc_overlap_matrix(OrbitalVector &braket);
ComplexMatrix calc_overlap_matrix(OrbitalVector &bra, OrbitalVector &ket);
ComplexMatrix calc_lowdin_matrix(OrbitalVector &Phi);

ComplexMatrix localize(double prec, OrbitalVector &Phi);
ComplexMatrix diagonalize(double prec, OrbitalVector &Phi, ComplexMatrix &F);
ComplexMatrix orthonormalize(double prec, OrbitalVector &Phi);

int size_empty(const OrbitalVector &vec);
int size_occupied(const OrbitalVector &vec);
int size_singly(const OrbitalVector &vec);
int size_doubly(const OrbitalVector &vec);
int size_paired(const OrbitalVector &vec);
int size_alpha(const OrbitalVector &vec);
int size_beta(const OrbitalVector &vec);
int get_multiplicity(const OrbitalVector &vec);
int get_electron_number(const OrbitalVector &vec, int spin = SPIN::Paired);

void set_spins(OrbitalVector &vec, const IntVector &spins);
void set_errors(OrbitalVector &vec, const DoubleVector &errors);
void set_occupancies(OrbitalVector &vec, const IntVector &occ);

IntVector get_spins(const OrbitalVector &vec);
IntVector get_occupancies(const OrbitalVector &vec);
DoubleVector get_errors(const OrbitalVector &vec);
DoubleVector get_norms(const OrbitalVector &vec);
DoubleVector get_squared_norms(const OrbitalVector &vec);

void print(const OrbitalVector &vec);

} //namespace orbital
} //namespace mrchem
