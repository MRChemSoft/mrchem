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

#include "DHScreening.h"
#include "GPESolver.h"
#include "Permittivity.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"

namespace mrchem {
/** @class PBESolver
 *  @brief class that performs the computation of the  ReactionPotential, named Self Consistent Reaction Field.
 */
class Nuclei;
class KAIN;
class PBESolver : public GPESolver {
public:
    PBESolver(Permittivity e,
              DHScreening kappa,
              const Nuclei &N,
              std::shared_ptr<mrcpp::PoissonOperator> P,
              std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
              double orb_prec,
              int kain_hist,
              int max_iter,
              bool acc_pot,
              bool dyn_thrs,
              std::string density_type);

    void setStaticSalt(bool do_static) { do_static_salt = do_static; }

    friend class ReactionPotential;

protected:
    bool do_static_salt;

    mrcpp::ComplexFunction kappa;
    mrcpp::ComplexFunction pbe_term;

    void setKappa(DHScreening k);

    // FIXME ComputeGamma should not be computing the PB term.
    //  THe PB term is actually an approximation for an external density, so
    //  it should be computed in the computedensities function or in the
    //  solvePoissonEquation function.
    void computeGamma(mrcpp::ComplexFunction &potential, mrcpp::ComplexFunction &out_gamma);
    void computePBTerm(mrcpp::ComplexFunction &V_tot, double salt_factor);
};
} // namespace mrchem
