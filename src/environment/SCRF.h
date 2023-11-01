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

#pragma once

#include "Permittivity.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"

namespace mrchem {
class KAIN;
/** @class SCRF
 *  @brief class that performs the computation of the  ReactionPotential, named Self Consistent Reaction Field.
 */
class SCRF final {
public:
    SCRF(const Permittivity &e,
         const Density &rho_nuc,
         std::shared_ptr<mrcpp::PoissonOperator> P,
         std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
         int kain_hist,
         int max_iter,
         bool dyn_thrs,
         const std::string &density_type);
    ~SCRF();

    double setConvergenceThreshold(double prec);

    mrcpp::ComplexFunction &getCurrentReactionPotential() { return this->Vr_n; }

    Permittivity &getPermittivity() { return this->epsilon; }

    void updateMOResidual(double const err_t) { this->mo_residual = err_t; }

    friend class ReactionPotential;
    auto computeEnergies(const Density &rho_el) -> std::tuple<double, double>;

protected:
    void clear();

private:
    bool dynamic_thrs;
    std::string density_type;

    int max_iter;
    int history;
    double apply_prec{-1.0};
    double conv_thrs{1.0};
    double mo_residual{1.0};

    Permittivity epsilon;

    Density rho_nuc; // As of right now, this is the biggest memory hog.
    // Alternative could be to precompute its contributions, as a potential is not as heavy as a density (maybe)
    // another one could be to define a representable function which only has the exact analytical form of the nuclear contribution.

    mrcpp::ComplexFunction Vr_n;

    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;
    std::shared_ptr<mrcpp::PoissonOperator> poisson;

    void computeDensities(OrbitalVector &Phi, Density &rho_out);
    void computeGamma(mrcpp::ComplexFunction &potential, mrcpp::ComplexFunction &out_gamma);

    mrcpp::ComplexFunction solvePoissonEquation(const mrcpp::ComplexFunction &ingamma, mrchem::OrbitalVector &Phi);

    void accelerateConvergence(mrcpp::ComplexFunction &dfunc, mrcpp::ComplexFunction &func, KAIN &kain);

    void nestedSCRF(const mrcpp::ComplexFunction &V_vac, std::shared_ptr<mrchem::OrbitalVector> Phi_p);
    mrcpp::ComplexFunction &setup(double prec, const std::shared_ptr<mrchem::OrbitalVector> &Phi_p);

    void resetComplexFunction(mrcpp::ComplexFunction &function);

    void printParameters() const;
    void printConvergenceRow(int i, double norm, double update, double time) const;
};
} // namespace mrchem
