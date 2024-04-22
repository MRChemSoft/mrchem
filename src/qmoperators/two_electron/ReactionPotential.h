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

#include "environment/GPESolver.h"
#include "qmoperators/QMPotential.h"

namespace mrchem {
/** @class ReactionPotential
 *  @brief class containing the solvent-substrate interaction reaction potential
 *  obtained by solving
 *  \f[
 *     \Delta V_{R} = -4\pi\left( \rho\frac{1-\epsilon}{\epsilon} + \gamma_s  \right)
 *  \f]
 *  where \f$\rho\f$ is the total molecular density of a solute molecule, \f$\epsilon\f$ is
 *  the Permittivity function of the continuum and \f$\gamma_s\f$ is the surface charge distribution.
 */
class ReactionPotential : public QMPotential {
public:
    /** @brief Initializes the ReactionPotential class.
     *  @param scrf A GPESolver instance which contains the parameters needed to compute the ReactionPotential.
     *  @param Phi A pointer to a vector which contains the orbitals optimized in the SCF procedure.
     */
    explicit ReactionPotential(std::unique_ptr<GPESolver> scrf, std::shared_ptr<mrchem::OrbitalVector> Phi = nullptr, bool mpi_share = false);
    ~ReactionPotential() override = default;

    GPESolver *getSolver() { return this->solver.get(); }

    /** @brief Updates the solver.mo_residual member variable. This variable is used to set the convergence criterion in
     * the dynamic convergence method. */
    void updateMOResidual(double const err_t) { this->solver->mo_residual = err_t; }

    friend class ReactionOperator;

protected:
    std::unique_ptr<GPESolver> solver;       //!< A GPESolver instance used to compute the ReactionPotential.
    std::shared_ptr<OrbitalVector> orbitals; ///< Unperturbed orbitals defining the ground-state electron density for the SCRF procedure.

    void setup(double prec) override;
    void clear() override;

private:
    virtual mrcpp::ComplexFunction &computePotential(double prec) const = 0;
};

} // namespace mrchem
