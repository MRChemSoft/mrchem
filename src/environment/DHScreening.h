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
#include "utils/print_utils.h"
#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>

namespace mrchem {
/** @class DHScreening
 *
 * @brief Square of the Debye-Huckel Screening parameter.
 * TODO write proper docs for this function
 *
 * Here the \f$\bar{\kappa}\f$ is dependent on a separate cavity function with its own radius.
 *
 */

class Cavity;

class DHScreening final : public mrcpp::RepresentableFunction<3> {
public:
    /** @brief Standard constructor. Initializes the #cavity_ion and #kappa_out with the input parameters.
     *  @param cavity_ion interlocking spheres of Cavity class.
     *  @param kappa_out value of the screening function outside the #cavity_ion.
     *  @param formulation Decides which formulation of the #DHScreening function to implement, only continuous screening function available
     * available as of now.
     */
    DHScreening(const Cavity cavity_ion, double kappa_out, std::string formulation);

    DHScreening(){};

    /** @brief Evaluates DHScreening at a point in 3D space with respect to the state of #inverse.
     *  @param r coordinates of a 3D point in space.
     *  @return \f$\frac{1}{\epsilon(\mathbf{r})}\f$ if #inverse is true, and \f$ \epsilon(\mathbf{r})\f$ if #inverse is
     *  false.
     */
    double evalf(const mrcpp::Coord<3> &r) const override;

    /** @brief Calls the Cavity::getCoordinates() method of the #cavity instance. */
    auto getCoordinates() const { return this->cavity_ion.getCoordinates(); }

    /** @brief Calls the Cavity::getRadii() method of the #cavity instance. */
    auto getRadii() const { return this->cavity_ion.getRadii(); }

    /** @brief Returns the value of #kappa_out. */
    auto getKOut() const { return this->kappa_out; }

    /** @brief Returns the cavity */
    Cavity getCavity() const { return this->cavity_ion; }

    /** @brief Returns the formulation */
    std::string getFormulation() const { return this->formulation; }

    /** @brief Print parameters */
    void printParameters() const;

private:
    double kappa_out;        //!< Dielectric constant describing the permittivity of the solvent.
    std::string formulation; //!< Formulation of the permittivity function, only exponential is used as of now.
    Cavity cavity_ion;       //!< A Cavity class instance.
};

} // namespace mrchem
