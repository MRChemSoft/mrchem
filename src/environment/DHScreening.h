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

#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>

#include "StepFunction.h"
#include "utils/print_utils.h"

namespace mrchem {
class Cavity;
/** @class DHScreening
 *
 * @brief Square of the Debye-Huckel Screening parameter.
 * @details
 * This is used for the Poisson-Boltzmann solver. The DHScreening function is defined as
 * \f[
 * \kappa^2(\mathbf{r}) = \begin{cases}
 * \kappa^2_{out} & \text{if } \mathbf{r} \notin \Omega_{ion} \\
 * 0.0 & \text{if } \mathbf{r} \in \Omega_{ion}
 * \end{cases}
 * \f]
 * This can be parametrized a number of ways. The one used here is
 * \f[
 * \kappa^2(\mathbf{r}) = (1 - C_{ion}(\mathbf{r})) \kappa^2_{out}
 * \f]
 * Where \f$C_{ion}(\mathbf{r})\f$ is the ion accessible Cavity function.
 */
class DHScreening final : public StepFunction {
public:
    /** @brief Standard constructor. Initializes the #cavity_ion and #kappa_out with the input parameters.
     *  @param cavity_ion interlocking spheres of Cavity class.
     *  @param kappa_out value of the screening function outside the #cavity_ion.
     *  @param formulation Decides which formulation of the #DHScreening function to implement, only continuous screening function available.
     * @details #kappa_out is given by
     * \f[
     * \kappa = \sqrt{\frac{2000 I_{0} e^2 N_a I_0}{\epsilon_{out} \epsilon_{in} k_B T}}
     * \f]
     * where \f$N_a\f$ is the Avogadro constant, e is the elementary charge, \f$I_0\f$ is the concentration of the ions,
     * \f$k_B\f$ is the Boltzmann constant, \f$T\f$ is the temperature, \f$\epsilon_{out}\f$ is the permittivity of the solvent and \f$\epsilon_{in}\f$ is the permittivity of free space.
     */
    DHScreening(std::shared_ptr<Cavity> cavity_ion, double kappa_out, const std::string &formulation);

    /** @brief Evaluates DHScreening at a point in 3D space.
     *  @param r coordinates of a 3D point in space.
     *  @return Value at point r.
     */
    double evalf(const mrcpp::Coord<3> &r) const override;
    void printParameters() const override { detail::print_header("Ion-accessible Cavity", this->formulation, getValueIn(), getValueOut()); }

private:
    std::string formulation{"Continuous Screening Function"}; //!< Formulation of the DHScreening function. Only linear variable is used now.
};

} // namespace mrchem
