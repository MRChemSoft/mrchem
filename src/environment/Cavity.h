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

#include <functional>
#include <vector>

#include <MRCPP/MWFunctions>
#include <MRCPP/mrcpp_declarations.h>

namespace mrchem {
/** @class Cavity
 * @brief Interlocking spheres cavity centered on the nuclei of the molecule.
 * The Cavity class represents the following function \cite Fosso-Tande2013
 *
 * \f[
 *    C(\mathbf{r})   = 1 - \prod^N_{i=1} (1-C_i(\mathbf{r})) \\
 *    C_i(\mathbf{r}) = 1 - \frac{1}{2}\left( 1 + \textrm{erf}\left(\frac{|\mathbf{r} - \mathbf{r}_i| -
 *    R_i}{\sigma_i}\right) \right)
 * \f]
 *
 * where \f$\mathbf{r}\f$ is the coordinate of a point in 3D space, \f$\mathbf{r}_i\f$
 * is the coordinate of the i-th nucleus, \f$R_i\f$ is the radius of the i-th sphere, and \f$\sigma_i\f$ is the width of
 * the transition between the inside and outside of the cavity. The transition has a sigmoidal shape, such that the
 * boundary is a smooth function instead of sharp boundaries often seen in other continuum models. This function is
 * \f$1\f$ inside and \f$0\f$ outside the cavity.
 *
 * The radii are computed as:
 *
 * \f[
 *   R_{i} =  \alpha_{i} R_{0,i} + \beta_{i}\sigma_{i}
 * \f]
 *
 * where:
 *
 *  - \f$R_{0,i}\f$ is the atomic radius. By default, the van der Waals radius.
 *  - \f$\alpha_{i}\f$ is a scaling factor. By default, 1.1
 *  - \f$\beta_{i}\f$ is a width scaling factor. By default, 0.5
 *  - \f$\sigma_{i}\f$ is the width. By default, 0.2
 */

class Cavity final : public mrcpp::RepresentableFunction<3> {
public:
    Cavity(const std::vector<mrcpp::Coord<3>> &coords, const std::vector<double> &R, const std::vector<double> &alphas, const std::vector<double> &betas, const std::vector<double> &sigmas);
    /** @brief Initializes the members of the class and constructs the analytical gradient vector of the Cavity.
     *
     * This CTOR applies a single width factor to the cavity and **does** not modify the radii. That is, in the formula:
     *
     * \f[
     *   R_{i} =  \alpha_{i} R_{0,i} + \beta_{i}\sigma_{i}
     * \f]
     *
     * for every atom \f$i\f$, \f$\alpha_{i} = 1.0\f$ and \f$\beta_{i} = 0.0\f$.
     */
    Cavity(const std::vector<mrcpp::Coord<3>> &coords, const std::vector<double> &R, double sigma)
            : Cavity(coords, R, std::vector<double>(R.size(), 1.0), std::vector<double>(R.size(), 0.0), std::vector<double>(R.size(), sigma)) {}
    double evalf(const mrcpp::Coord<3> &r) const override;
    auto getGradVector() const { return this->gradvector; }

    /** @brief Print parameters */
    void printParameters() const;

protected:
    std::vector<double> radii_0;                                             //!< Contains the *unscaled* radius of each sphere in #Center.
    std::vector<double> alphas;                                              //!< The radius scaling factor for each sphere.
    std::vector<double> betas;                                               //!< The width scaling factor for each sphere.
    std::vector<double> sigmas;                                              //!< The width for each sphere.
    std::vector<double> radii;                                               //!< Contains the radius of each sphere in #Center. \f$R_i = \alpha_{i} R_{0,i} + \beta_{i}\sigma_{i}\f$
    std::vector<mrcpp::Coord<3>> centers;                                    //!< Contains each of the spheres centered on the nuclei of the Molecule.
    std::vector<std::function<double(const mrcpp::Coord<3> &r)>> gradvector; //< Analytical derivatives of the Cavity.

    bool isVisibleAtScale(int scale, int nQuadPts) const override;
    bool isZeroOnInterval(const double *a, const double *b) const override;
};

namespace detail {
auto gradCavity(const mrcpp::Coord<3> &r, int index, const std::vector<mrcpp::Coord<3>> &centers, const std::vector<double> &radii, const std::vector<double> &width) -> double;
} // namespace detail
} // namespace mrchem
