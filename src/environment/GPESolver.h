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

#include <memory>
#include <tuple>

#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>

#include "DHScreening.h"
#include "Permittivity.h"
#include "qmfunctions/Density.h"

namespace mrchem {
class KAIN;
class ReactionPotential;
class ReactionPotentialD1;
class ReactionPotentialD2;

enum class SCRFDensityType : int { TOTAL = 0, ELECTRONIC = 1, NUCLEAR = 2 };

/** @class GPESolver
 *  @brief Solves the Generalized Poisson equation iteratively
 *  @details The Generalized Poisson equation is given by
 * \f[
 * \nabla \cdot \left( \epsilon(\mathbf{r}) \nabla V(\mathbf{r}) \right) = -4\pi \rho(\mathbf{r})
 * \f]
 * where \f$\epsilon(\mathbf{r})\f$ is the permittivity, \f$V(\mathbf{r})\f$ is the total electrostatic potential and \f$\rho(\mathbf{r})\f$ is the molecular charge density defined as:
 * \f[
 * \rho(\mathbf{r}) = \rho_{el}(\mathbf{r}) + \rho_{nuc}(\mathbf{r})
 * \f]
 * where \f$\rho_{el}\f$ is the electronic charge density and \f$\rho_{nuc}\f$ is the nuclear charge density.
 * The Generalized Poisson equation is solved iteratively through a set of micro-iteration on each SCF-iteration by appliation of the Poisson operator \f$\mathcal{P}\f$  :cite:`Fosso-Tande2013`
 * \f[
 * V_R(\mathbf{r}) = \mathcal{P} \star \left[ \rho_{eff}(\mathbf{r}) - \rho(\mathbf{r}) + \gamma_s(\mathbf{r}) \right]
 * \f]
 * where \f$\gamma_s(\mathbf{r})\f$ is the surface charge distribution describing the polarization at the surface, \f$\rho_{eff}(\mathbf{r})\f$ is the effective charge density given by
 * \f$\frac{\rho(\mathbf{r})}{\epsilon(\mathbf{r})}\f$ and \f$V_R(\mathbf{r})\f$ is the reaction potential.
 *
 * We utilize a so-called dynamic threshold to more easily converge the reaction potential. This is done by setting the convergence threshold of the micro-iterations to
 * the MO update of the previous SCF iteration, unless the MO update is small enough (once the quality of the MOs is good enough, we use the default convergence threshold).
 * Another optimization used is that we utilize the previous SCF converged Reaction potential as an initial guess for the next micro-iterations. These procedures are
 * investigated and explained in :cite:`gerez2023`
 */
class GPESolver {
public:
    GPESolver(const Permittivity &e,
              const Density &rho_nuc,
              std::shared_ptr<mrcpp::PoissonOperator> P,
              std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
              int kain_hist,
              int max_iter,
              bool dyn_thrs,
              SCRFDensityType density_type);
    ~GPESolver();

    /** @brief Sets the convergence threshold for the micro-iterations, used with dynamic thresholding.
     *  @param prec value to set the convergence threshold to
     *  @return the current convergence threshold.
     *  @details will check if the MO update is small enough (ten times as big) wrt. to the scf convergence threshold, if so, it will use the default convergence threshold.
     * If not, it will use the MO update as the convergence threshold.
     */
    double setConvergenceThreshold(double prec);

    Permittivity &getPermittivity() { return this->epsilon; }

    void updateMOResidual(double const err_t) { this->mo_residual = err_t; }

    /** @brief Computes the energy contributions from the reaction potential
     * @param rho_el the electronic charge density
     * @return a tuple containing the electronic and nuclear energy contributions
     * @details We compute the reaction energy through the following integral:
     * \f[
     * E_{R} = \frac{1}{2}\int \rho_{el}(\mathbf{r}) V_R(\mathbf{r}) d\mathbf{r} + \frac{1}{2} \int \rho_{nuc}(\mathbf{r}) V_R(\mathbf{r}) d\mathbf{r}
     * \f]
     * Each term represents the electronic and nuclear contributions to the reaction energy, respectively.
     * We compute each term separately, and return a tuple containing both.
     */
    auto computeEnergies(const Density &rho_el) -> std::tuple<double, double>;

    auto getDensityType() const -> SCRFDensityType { return this->density_type; }

    friend class ReactionPotential;
    friend class ReactionPotentialD1;
    friend class ReactionPotentialD2;

protected:
    bool dynamic_thrs;
    SCRFDensityType density_type; //!< Decides which density we will use for computing the reaction potential, options are ``total``, ``electronic`` and ``nuclear``.

    int max_iter;
    int history;
    double apply_prec{-1.0};
    double conv_thrs{1.0};
    double mo_residual{1.0};
    std::string solver_name{"Generalized Poisson"};

    Permittivity epsilon;

    Density rho_nuc; // As of right now, this is the biggest memory hog.
    // Alternative could be to precompute its contributions, as a potential is not as heavy as a density (maybe)
    // another one could be to define a representable function which only has the exact analytical form of the nuclear contribution.
    // Since we already have \nabla^2 V_nuc = -4*pi*rho_nuc (nuclear potential) we could use this as a way to bypass computing rho_nuc at all
    // Same with the coulomb potential, which is basically what is left from V_vac after subtracting V_nuc. in one way we could just precompute both and
    // just iterate through V_R only. Only issue here is (1 -1\epsilon)/\epsilon * \rho_nuc as I am not sure how to represent this as an analytitcal function,
    // maybe after applying the Poisson operator?

    mrcpp::ComplexFunction Vr_n;

    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;
    std::shared_ptr<mrcpp::PoissonOperator> poisson;

    void clear();

    /** @brief computes density wrt. the density_type variable
     * @param Phi the molecular orbitals
     * @param rho_out Density function in which the density will be computed.
     * @details The total charge density is given by the sum of the electronic and nuclear charge densities:
     * \f[
     * \rho(\mathbf{r}) = \rho_{el}(\mathbf{r}) + \rho_{nuc}(\mathbf{r})
     * \f]
     * where \f$\rho_{el}\f$ is the electronic charge density and \f$\rho_{nuc}\f$ is the nuclear charge density.
     * The nuclear charge density is stored in the class variable #rho_nuc, while we compute the electronic charge density
     * from the molecular orbitals.
     * The class variable #density_type decides the density which will be computed in rho_out, options are ``total``, ``electronic`` and ``nuclear``.
     */
    void computeDensities(const Density &rho_el, Density &rho_out);

    /** @brief Computes the surface charge distibution due to polarization at the solute-solvent boundary
     * @param potential Potential used to compute \f$\nabla V(\mathbf{r})\f$
     * @param out_gamma ComplexFunction in which the surface charge distribution will be computed.
     * @details The surface charge distribution is given by
     * \f[
     * \gamma_s(\mathbf{r}) = \frac{\log \frac{\epsilon_{in}}{\epsilon_{out}}}{4 \pi } \left( \nabla C(\mathbf{r}) \cdot \nabla V(\mathbf{r})\right)
     * \f]
     * where \f$\epsilon_{in}\f$ is the permittivity inside the cavity and \f$\epsilon_{out}\f$ is the permittivity outside the cavity.
     */
    virtual void computeGamma(mrcpp::ComplexFunction &potential, mrcpp::ComplexFunction &out_gamma);

    /** @brief Iterates once through the Generalized Poisson equation to compute the reaction potential
     * @param ingamma the surface charge distribution
     * @param Phi the molecular orbitals
     * @return the reaction potential
     * @details Constructs the effective charge density \f$\rho_{eff}(\mathbf{r})\f$ and the Poisson operator \f$\mathcal{P}\f$ as:
     * \f[
     * V_R(\mathbf{r}) = \mathcal{P} \star \left[ \rho_{eff}(\mathbf{r}) - \rho(\mathbf{r}) + \gamma_s(\mathbf{r}) \right]
     * \f]
     * where \f$\gamma_s(\mathbf{r})\f$ is the surface charge distribution describing the polarization at the surface, \f$\rho_{eff}(\mathbf{r})\f$ is the effective charge density given by
     * \f$\frac{\rho(\mathbf{r})}{\epsilon(\mathbf{r})}\f$ and \f$V_R(\mathbf{r})\f$ is the reaction potential.
     */
    mrcpp::ComplexFunction solvePoissonEquation(const mrcpp::ComplexFunction &ingamma, const Density &rho_el);

    /** @brief Uses KAIN to accelerate convergece of the reaction potential
     * @param dfunc the current update of the reaction potential
     * @param func the current reaction potential
     * @param kain the KAIN object
     */
    void accelerateConvergence(mrcpp::ComplexFunction &dfunc, mrcpp::ComplexFunction &func, KAIN &kain);

    /** @brief Iterates through the application of the Poisson operator to Solve the Generalized Poisson equation
     *  @param V_vac the vacuum potential
     *  @param Phi_p the molecular orbitals
     * @details Iterating through the application of the Poisson operator is done through a set of micro-iterations,
     * where the convergence threshold is set to the MO update of the previous SCF iteration.
     * The micro-iterations are done through the following steps:
     *  -# Compute the total potential as \f$V(\mathbf{r}) = V_{vac}(\mathbf{r}) + V_R(\mathbf{r})\f$
     *  -# Compute the surface charge distribution \f$\gamma_s(\mathbf{r})\f$ with #computeGamma
     *  -# Compute a new reaction potential \f$V_R(\mathbf{r})\f$ by application of the Poisson operator with #solvePoissonEquation
     *  -# Calculate the update of the reaction potential as \f$\Delta V_R(\mathbf{r}) = V_R(\mathbf{r}) - V_R^{old}(\mathbf{r})\f$
     *  -# Accelerate convergence of the reaction potential through KAIN
     *  -# Update the reaction potential as \f$V_R(\mathbf{r}) = V_R^{old}(\mathbf{r}) + \Delta V_R(\mathbf{r})\f$
     *  -# Check if the reaction potential has converged, if not, repeat from step 1.
     */
    void runMicroIterations(const mrcpp::ComplexFunction &V_vac, const Density &rho_el);

    /** @brief Setups and computes the reaction potential through the microiterations
     * @param V_vac the vacuum potential
     * @param Phi_p the molecular orbitals
     * @return The converged reaction potential for the current SCF iteration
     * @details An initial guess of the reaction potential is computed with the following steps:
     * -# Set the total potential as \f$V(\mathbf{r}) = V_{vac}(\mathbf{r})\f$
     * -# Compute the surface charge distribution \f$\gamma_s(\mathbf{r})\f$ from this potential
     * -# Iterate once through the application of the Poisson operator to return the initial guess of the reaction potential \f$V_R(\mathbf{r})\f$
     *
     * the method then runs the micro-iterations through #runMicroIterations and returns the converged reaction potential.
     * If this is not the first SCF iteration, the previous converged reaction potential is used as an initial guess for the micro-iterations.
     */
    mrcpp::ComplexFunction &solveEquation(double prec, const Density &rho_el);

    /** @brief Frees the memory used by the FunctionTrees of the input Complexfunction and reallocates them.
     * @param function the ComplexFunction to reset
     * @details This is done to avoid memory leaks when the ComplexFunction is used in the micro-iterations.
     */
    void resetComplexFunction(mrcpp::ComplexFunction &function);

    virtual void printParameters() const;
    void printConvergenceRow(int i, double norm, double update, double time) const;
};
} // namespace mrchem
