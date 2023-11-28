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
 *  @brief Solves the Poisson-Boltzmann equation iteratively
 *  @details The Poisson-Boltzmann equation is solved iteratively using the SCRF procedure outlined in GPESolver.
 * The Poisson-Boltzmann equation models the electrostatic potential in a solvent with electrolytes.
 * The general equation for electrolyte solutions is given by
 * \f[
 * \nabla \cdot \epsilon \nabla V_{tot} = -4\pi\left(\rho_{el} + \rho_{nuc} + \rho_{ext}\right)
 * \f]
 * where \f$ V_{tot} \f$ is the total electrostatic potential, \f$ \epsilon\f$ is the permittivity function of the solvent,
 * \f$\rho_{el}\f$ is the electronic charge density, \f$\rho_{nuc}\f$ is the nuclear charge density and \f$\rho_{ext}\f$ is the external charge density.
 * In the general form for the Poisson-Boltzmann equation, the external charge density is approximated by assuming a boltzmann distribution of the ions.
 * \f[
 * \rho_{ext} = \sum_{i}^{N_{ion}} q_i e I_{0, i}\exp\left(-\frac{q_i e V_{tot}}{k_B T}\right)
 * \f]
 * where \f$I_{0, i}\f$ is the concentration of the i-th ion species, \f$q_i\f$ is the charge of the i-th ion species, \f$k_B\f$ is the Boltzmann constant and \f$T\f$ is the temperature.
 * In this implementation we assume a 1:1 (\f$I_{0, 0} = I_{0, 1}\f$) electrolyte soluttion of ions of same opposite charges (\f$z_i = +1, -1\f$). This simplifies the external density to
 * \f[
 * \rho_{ext} = -2 e I_{0} \sinh\left(\frac{e V_{tot}}{2 k_B T}\right)
 * \f]
 * where \f$I_{0}\f$ is the concentration of the ions.
 * We can plug this into the first equation (and massage terms a bit) to arrive at the Poisson-Boltzmann equation for 1:1 electrolyte solution
 * \f[
 * \nabla^2 V_{R} = -4\pi\frac{1-\epsilon}{\epsilon}\left(\rho_{el} + \rho_{nuc}\right) + \gamma_s - \kappa^2 \sinh\left(V_{tot}\right)
 * \f]
 * where \f$\gamma_s\f$ is the surface charge density, \f$\kappa\f$ is obtained from the DHScreening class and \f$V_{R}\f$ is the reaction potential.
 */
class PBESolver : public GPESolver {
public:
    PBESolver(const Permittivity &e,
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
    DHScreening kappa; ///< the DHScreening object used to compute the PB term \f$\kappa\f$
    std::string solver_name{"Poisson-Boltzmann"};

    /** @brief constructs the surface chage distribution and adds it to the PB term
     * @param[in] potential the potential to compute \f$\nabla V\f$ from
     * @param[out] out_gamma the ComplexFunction in which to store the result
     * @details Method follows the implementation in GPESolver::computeGamma, but adds the PB term to the surface charge distribution.
     */
    void computeGamma(mrcpp::ComplexFunction &potential, mrcpp::ComplexFunction &out_gamma) override;

    /** @brief Computes the PB term
     * @param[in] V_tot the total potential
     * @param[in] salt_factor the salt factor deciding how much of the total concentration to include in the PB term
     * @param[out] pb_term the ComplexFunction in which to store the result
     * @details The PB term is computed as \f$ \kappa^2 \sinh(V_{tot}) \f$ and returned.
     */
    virtual void computePBTerm(mrcpp::ComplexFunction &V_tot, const double salt_factor, mrcpp::ComplexFunction &pb_term);
};
} // namespace mrchem
