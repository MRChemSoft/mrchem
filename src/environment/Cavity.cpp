/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
#include "Cavity.h"
#include "chemistry/PhysicalConstants.h"
#include "utils/math_utils.h"

namespace mrchem {

/** @brief Initializes the members of the class and constructs the analytical gradient vector of the Cavity. */
Cavity::Cavity(std::vector<mrcpp::Coord<3>> &centers, std::vector<double> &radii, double width)
        : width(width)
        , radii(radii)
        , centers(centers) {
    setGradVector();
}

/**  @relates mrchem::Cavity
 *   @brief Constructs a single element of the gradient of the Cavity.
 *
 *   This constructs the analytical partial derivative of the Cavity \f$C\f$ with respect to \f$x\f$, \f$y\f$ or \f$z\f$
 * coordinates and evaluates it at a point \f$\mathbf{r}\f$. This is given for \f$x\f$ by
 * \f[
 *    \frac{\partial C\left(\mathbf{r}\right)}{\partial x} = \left(1 - C{\left(\mathbf{r} \right)}\right)
 *                                                           \sum_{i=1}^{N} - \frac{\left(x-{x}_{i}\right)e^{-
 * \frac{\operatorname{s_{i}}^{2}{\left(\mathbf{r} \right)}}{\sigma^{2}}}}
 *                                                           {\sqrt{\pi}\sigma\left(0.5
 * \operatorname{erf}{\left(\frac{\operatorname{s_{i}}{\left(\mathbf{r}  \right)}}{\sigma} \right)}
 *                                                           + 0.5\right) \left| \mathbf{r} - \mathbf{r}_{i} \right|}
 * \f]
 * where the subscript \f$i\f$ is the index related to each
 * sphere in the cavity, and \f$\operatorname{s}\f$ is the signed normal distance from the surface of each sphere.
 *   @param r The coordinates of a test point in 3D space.
 *   @param index An integer that defines the variable of differentiation (0->x, 1->z and 2->z).
 *   @param centers A vector containing the coordinates of the centers of the spheres in the cavity.
 *   @param radii A vector containing the radii of the spheres.
 *   @param width A double value describing the width of the transition at the boundary of the spheres.
 *   @return A double number which represents the value of the differential (w.r.t. x, y or z) at point r.
 */
auto gradCavity(const mrcpp::Coord<3> &r, int index, const std::vector<mrcpp::Coord<3>> &centers, std::vector<double> &radii, double width) -> double {
    double C = 1.0;
    double DC = 0.0;
    double sqrt_pi = std::sqrt(mrcpp::pi);
    for (int i = 0; i < centers.size(); i++) {
        double s = math_utils::calc_distance(centers[i], r) - radii[i];
        double ds = (r[index] - centers[i][index]) / (math_utils::calc_distance(centers[i], r));
        double Theta = 0.5 * (1 + std::erf(s / width));
        double Ci = 1.0 - Theta;
        C *= 1.0 - Ci;

        double DCi = -(1.0 / (width * sqrt_pi)) * std::exp(-std::pow(s / width, 2.0)) * ds;

        double numerator = DCi;
        double denominator = 1.0 - Ci;

        if (((1.0 - Ci) < 1.0e-12) and ((1.0 - Ci) >= 0)) {
            denominator = 1.0e-12;
        } else if (((1.0 - Ci) > -1.0e-12) and ((1.0 - Ci) <= 0)) {
            denominator = -1.0e-12;
        }

        if ((DCi < 1.0e-12) and (DCi >= 0)) {
            numerator = 1.0e-12;
        } else if ((DCi > -1.0e-12) and (DCi <= 0)) {
            numerator = -1.0e-12;
        }
        DC += numerator / denominator;
    }
    DC = C * DC;
    return DC;
}
/** @brief Sets the different partial derivatives in the \link #gradvector gradient \endlink of the Cavity. */
void Cavity::setGradVector() {
    auto p_gradcavity = [this](const mrcpp::Coord<3> &r, int index) { return gradCavity(r, index, centers, radii, width); };
    for (auto i = 0; i < 3; i++) {
        this->gradvector.push_back([i, p_gradcavity](const mrcpp::Coord<3> &r) -> double { return p_gradcavity(r, i); });
    }
}
/** @brief Evaluates the value of the cavity at a 3D point \f$\mathbf{r}\f$
 *  @param r coordinate of 3D point at which the Cavity is to be evaluated at.
 *  @return double value of the Cavity at point \f$\mathbf{r}\f$
 */
double Cavity::evalf(const mrcpp::Coord<3> &r) const {
    double C = 1.0;
    for (int i = 0; i < centers.size(); i++) {
        double s = math_utils::calc_distance(centers[i], r) - radii[i];
        double Theta = 0.5 * (1 + std::erf(s / width));
        double Ci = 1 - Theta;
        C *= 1 - Ci;
    }
    C = 1 - C;
    return C;
}

bool Cavity::isVisibleAtScale(int scale, int nQuadPts) const {

    auto visibleScale = static_cast<int>(-std::floor(std::log2(nQuadPts * 2.0 * this->width)));

    if (scale < visibleScale) { return false; }

    return true;
}

bool Cavity::isZeroOnInterval(const double *a, const double *b) const {
    for (int k = 0; k < centers.size(); k++) {
        for (int i = 0; i < 3; i++) {
            double cavityMinOut = (this->centers[k][i] - radii[i]) - 3.0 * this->width;
            double cavityMinIn = (this->centers[k][i] - radii[i]) + 3.0 * this->width;
            double cavityMaxIn = (this->centers[k][i] + radii[i]) - 3.0 * this->width;
            double cavityMaxOut = (this->centers[k][i] + radii[i]) + 3.0 * this->width;
            if (a[i] > cavityMaxOut or (a[i] > cavityMinIn and b[i] < cavityMaxIn) or b[i] < cavityMinOut) { return true; }
        }
    }
    return false;
}

} // namespace mrchem
