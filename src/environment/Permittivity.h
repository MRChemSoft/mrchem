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

#include "Cavity.h"
#include "StepFunction.h"
#include "utils/print_utils.h"
#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>

namespace mrchem {
/** @class Permittivity
 *
 * @brief Permittivity function related to a substrate molecule and a solvent continuum.
 * The Permittivity class represents the following function \cite Fosso-Tande2013
 * \f[
 *   \epsilon(\mathbf{r}) = \epsilon_{in}\exp\left(\left(\log\frac{\epsilon_{out}}{\epsilon_{in}} \right)
 *                                               \left(1 - C(\mathbf{r})\right)\right)
 * \f]
 * where \f$\mathbf{r}\f$ is the coordinate of a point in 3D space, \f$ C \f$ is the #cavity function of the substrate,
 * and \f$\epsilon_{in}\f$ and \f$ \epsilon_{out} \f$ are the dielectric constants describing, respectively, the
 * permittivity \link #epsilon_in inside \endlink  and \link #epsilon_out outside\endlink the #cavity of the substrate.
 */

class Cavity;

class Permittivity final : public StepFunction {
public:
    /** @brief Standard constructor. Initializes the #cavity, #epsilon_in and #epsilon_out with the input parameters.
     *  @param cavity interlocking spheres of Cavity class.
     *  @param epsilon_in permittivity inside the #cavity.
     *  @param epsilon_out permittivity outside the #cavity.
     *  @param formulation Decides which formulation of the #Permittivity function to implement, only exponential
     * available as of now.
     */
    Permittivity(std::shared_ptr<Cavity> cavity, double epsilon_in, double epsilon_out, std::string formulation);

    /** @brief Evaluates Permittivity at a point in 3D space with respect to the state of #inverse.
     *  @param r coordinates of a 3D point in space.
     *  @return \f$\frac{1}{\epsilon(\mathbf{r})}\f$ if #inverse is true, and \f$ \epsilon(\mathbf{r})\f$ if #inverse is
     *  false.
     */
    double evalf(const mrcpp::Coord<3> &r) const override;

    void printHeader() const override { detail::print_header("Solvation Cavity", this->formulation, getValueIn(), getValueOut()); }

private:
    std::string formulation{"exponential"};
};

} // namespace mrchem
