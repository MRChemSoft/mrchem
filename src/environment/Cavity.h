/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <MRCPP/MWFunctions>

namespace mrchem {
/** @class Cavity
 * @brief Interlocking spheres cavity centered on the nuclei of the molecule.
 * The Cavity class represents the following function \cite Fosso-Tande2013
 * \f[
 *    C(\mathbf{r})   = 1 - \prod^N_{i=1} (1-C_i(\mathbf{r})) \\
 *    C_i(\mathbf{r}) = 1 - \frac{1}{2}\left( 1 + \textrm{erf}\left(\frac{|\mathbf{r} - \mathbf{r}_i| -
 *    R_i}{\sigma}\right) \right)
 * \f]
 * where \f$\mathbf{r}\f$ is the coordinate of a point in 3D space, \f$\mathbf{r}_i\f$
 * is the coordinate of the i-th nucleus, \f$R_i\f$ is the radius of the i-th sphere, and \f$\sigma\f$ is the #width of
 * the transition between the inside and outside of the cavity. The transition has a sigmoidal shape, such that the
 * boundary is a smooth function instead of sharp boundaries often seen in other continuum models. This function is
 * \f$1\f$ inside and \f$0\f$ outside the cavity.
 */

class Cavity final : public mrcpp::RepresentableFunction<3> {
public:
    Cavity(std::vector<mrcpp::Coord<3>> &centers, std::vector<double> &radii, double width);
    double evalf(const mrcpp::Coord<3> &r) const override;
    auto getGradVector() const { return this->gradvector; }
    std::vector<mrcpp::Coord<3>> getCoordinates() const { return centers; } //!< Returns #centers.
    std::vector<double> getRadii() const { return radii; }                  //!< Returns #radii.
    friend class Permittivity;

protected:
    double width;                         //!< width of the Cavity boundary.
    std::vector<double> radii;            //!< Contains the radius of each sphere in #Center.
    std::vector<mrcpp::Coord<3>> centers; //!< Contains each of the spheres centered on the nuclei of the Molecule.
    std::vector<std::function<double(const mrcpp::Coord<3> &r)>> gradvector; //< Analytical derivatives of the Cavity.

    void setGradVector();
    bool isVisibleAtScale(int scale, int nQuadPts) const override;
    bool isZeroOnInterval(const double *a, const double *b) const override;
};
} // namespace mrchem
