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
 *  @brief Solves the Linearized Poisson-Boltzmann equation iteratively
 * @details The Linearized Poisson-Boltzmann equation is solved iteratively using the SCRF procedure outlined in GPESolver and PBESolver.
 * The linearized Poisson-Boltzmann equation is a further simplification of the Poisson-Boltzmann equation, outlined in PBESolver, where the PB term is expanded and
 * only the linear term is included. This is a good approximation for low ionic strength solutions.
 * The linearized Poisson-Boltzmann equation is given by
 * \f[
 * \nabla^2 V_{R} = -4\pi\frac{1-\epsilon}{\epsilon}\left(\rho_{el} + \rho_{nuc}\right) + \gamma_s - \kappa^2 V_{tot}
 * \f]
 * where \f$\gamma_s\f$ is the surface charge density, \f$\kappa\f$ is obtained from the DHScreening class and \f$V_{R}\f$ is the reaction potential.
 */
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
               SCRFDensityType density_type);

    friend class ReactionPotential;

protected:
    /** @brief Computes the PB term
     * @param[in] V_tot the total potential
     * @param[in] salt_factor the salt factor deciding how much of the total concentration to include in the PB term
     * @param[out] pb_term the ComplexFunction in which to store the result
     * @details The PB term is computed as \f$ \kappa^2 V_{tot} \f$ and returned.
     */
    void computePBTerm(mrcpp::ComplexFunction &V_tot, const double salt_factor, mrcpp::ComplexFunction &pb_term) override;
    std::string solver_name{"Linearized Poisson-Boltzmann"};
};
} // namespace mrchem
