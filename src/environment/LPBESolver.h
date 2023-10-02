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
#include "PBESolver.h"
#include "Permittivity.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"

namespace mrchem {
/** @class LPBESolver
 *  @brief class that performs the computation of the  ReactionPotential, named Self Consistent Reaction Field.
 */
class Nuclei;
class KAIN;
class LPBESolver final : public PBESolver {
public:
    LPBESolver(const Permittivity &e,
               const DHScreening &k,
               const Density &rho_nuc,
               std::shared_ptr<mrcpp::PoissonOperator> P,
               std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
               int kain_hist,
               int max_iter,
               bool dyn_thrs,
               const std::string &density_type)
            : PBESolver(e, kappa, N, P, D, orb_prec, kain_hist, max_iter, acc_pot, dyn_thrs, density_type) {}

    friend class ReactionPotential;

protected:
    void computePBTerm(mrcpp::ComplexFunction &V_tot, const double salt_factor, mrcpp::ComplexFunction &pb_term) override;
};
} // namespace mrchem
