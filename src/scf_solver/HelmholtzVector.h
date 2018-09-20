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

#include "qmfunctions/qmfunction_fwd.h"
#include "qmoperators/qmoperator_fwd.h"

/** @class HelmholtzVector
 *
 * @brief Container of HelmholtzOperators for a corresponding OrbtialVector
 *
 * This class assigns one HelmholtzOperator to each orbital in an OrbitalVector.
 * The operators can be re-used if several orbitals share the same (or similar)
 * energy, or if the change in energy is small relative to the previous iteration.
 */

namespace mrchem {

class HelmholtzVector final {
public:
    HelmholtzVector(double build, double thrs = -1.0);

    void setup(double prec, const DoubleVector &energies);
    void clear();

    void setThreshold(double thrs) { this->threshold = thrs; }
    double getThreshold() const { return this->threshold; }

    double getLambda(int i) const { return this->lambda[i]; }
    DoubleVector getLambdaVector() const;
    ComplexMatrix getLambdaMatrix() const;

    mrcpp::HelmholtzOperator &operator[](int i);
    const mrcpp::HelmholtzOperator &operator[](int i) const;

    int printTreeSizes() const;

    Orbital operator()(int i, Orbital inp);
    OrbitalVector operator()(OrbitalVector &inp);

private:
    double threshold;  ///< For re-using operators. Negative means always recreate
    double build_prec; ///< Precision for construction of Helmholtz operators
    double apply_prec; ///< Precision for application of Helmholtz operators

    std::vector<int> oper_idx;  ///< Points to a HelmholtzOperator in the operators vector
    std::vector<double> lambda; ///< The lambda value used for the corresponding HelmholtzOperator
    std::vector<mrcpp::HelmholtzOperator *> operators; ///< Vector of Helmholtz operators

    int initHelmholtzOperator(double energy, int i);
    void clearUnused();
};

} // namespace mrchem
