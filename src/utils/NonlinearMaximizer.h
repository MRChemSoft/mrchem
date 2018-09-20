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

#include <Eigen/Core>

#include "mrchem.h"

/* @class NonlinearMaximizer
 *
 * @date Jan 31, 2013
 * @author Peter Wind <peter.wind@uit.no> \n
 *         CTCC, University of Tromsø
 *
 * @brief Maximization of nonlinear functional
 */

namespace mrchem {

class NonlinearMaximizer {
public:
    NonlinearMaximizer() : N2h(0) { }
    int maximize();

protected:
    int N2h;//size (for orbital localization: N2h = N*(N-1)/2)
    DoubleMatrix hessian;
    DoubleVector gradient;

    virtual double functional() const { return 0.0; }
    virtual double make_gradient() { return -1.0; }
    virtual double make_hessian() { return -1.0; }
    virtual double get_hessian(int i, int j) { return -1.0; }
    virtual void multiply_hessian(DoubleVector &vec, DoubleVector &Hv){}
    virtual void do_step(const DoubleVector &step) { }
};

} //namespace mrchem
