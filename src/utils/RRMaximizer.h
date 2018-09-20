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

#include "utils/NonlinearMaximizer.h"
#include "qmfunctions/qmfunction_fwd.h"

/** subclass which defines the particular Gradient and Hessian
 * and other specific functions for a maximization of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 *
 */

namespace mrchem {

class RRMaximizer final : public NonlinearMaximizer {
public:
    RRMaximizer(double prec, OrbitalVector &Phi);
    const DoubleMatrix &getTotalU() const { return this->total_U; }
    double get_hessian(int i, int j);
    void multiply_hessian(DoubleVector &vec, DoubleVector &Hv);

protected:
    int N;                  // number of orbitals
    DoubleMatrix r_i_orig;  // <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
    DoubleMatrix r_i ;      // rotated  r_i_orig
    DoubleMatrix total_U;   // the rotation matrix of the orbitals

    double functional() const;
    double make_gradient();
    double make_hessian();
    void do_step(const DoubleVector &step);
};

} //namespace mrchem
