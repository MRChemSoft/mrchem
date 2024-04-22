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

#include "Cavity.h"
#include "utils/print_utils.h"

namespace mrchem {

namespace detail {
void print_header(const std::string &header, const std::string &formulation, double in_value, double out_value);
} // namespace detail

class Cavity;

/** @class StepFunction
 *
 * @brief StepFunction function related to a Cavity function.
 * The StepFunction class represents the following equation
 * \f[
 *   S(\mathbf{r}) = \begin{cases} S_{in} & \text{if } \mathbf{r} inside C \\
 *                     S_{out} & \text{if } \mathbf{r} outside C
 * \end{cases}
 * \f]
 * where \f$\mathbf{r}\f$ is the coordinate of a point in 3D space, \f$ C \f$ is the #cavity function,
 * and \f$S_{in}\f$ and \f$ S_{out} \f$ are the values of the function inside and outside the #cavity respectively.
 */
class StepFunction : public mrcpp::RepresentableFunction<3> {
public:
    /** @brief Standard constructor. Initializes the #cavity, #in and #out with the input parameters.
     *  @param cavity interlocking spheres of Cavity class.
     *  @param epsilon_in permittivity inside the #cavity.
     *  @param epsilon_out permittivity outside the #cavity.
     * available as of now.
     */
    StepFunction(std::shared_ptr<Cavity> cavity, double val_in, double val_out);

    auto getValueIn() const { return this->in; }
    auto getValueOut() const { return this->out; }

    std::shared_ptr<Cavity> getCavity_p() const { return this->cavity_p; }

    virtual void printParameters() const { detail::print_header("Step function", "Standard", getValueIn(), getValueOut()); }

protected:
    double in;                        //!< Value of the function inside the #cavity.
    double out;                       //!< value of the function outside the #cavity.
    std::shared_ptr<Cavity> cavity_p; //!< A shared pointer to a Cavity class instance.
};

} // namespace mrchem
