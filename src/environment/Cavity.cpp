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
#include "Cavity.h"

#include <cmath>
#include <functional>
#include <vector>

#include <MRCPP/mrcpp_declarations.h>

#include "chemistry/PhysicalConstants.h"
#include "utils/math_utils.h"

namespace mrchem {
namespace detail {
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
auto gradCavity(const mrcpp::Coord<3> &r, int index, const std::vector<mrcpp::Coord<3>> &centers, const std::vector<double> &radii, const std::vector<double> &widths) -> double {
    auto C = 1.0;
    auto DC = 0.0;
    auto sqrt_pi = std::sqrt(mrcpp::pi);

    for (int i = 0; i < centers.size(); ++i) {
        auto center = centers[i];
        auto radius = radii[i];
        auto sigma = widths[i];

        auto s = math_utils::calc_distance(center, r) - radius;
        auto ds = (r[index] - center[index]) / (math_utils::calc_distance(center, r));
        auto Theta = 0.5 * (1 + std::erf(s / sigma));
        auto Ci = 1.0 - Theta;
        C *= 1.0 - Ci;

        double DCi = -(1.0 / (sigma * sqrt_pi)) * std::exp(-std::pow(s / sigma, 2.0)) * ds;

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
} // namespace detail

/** @brief Initializes the members of the class and constructs the analytical gradient vector of the Cavity. */
Cavity::Cavity(const std::vector<mrcpp::Coord<3>> &coords, const std::vector<double> &R, const std::vector<double> &A, const std::vector<double> &B, const std::vector<double> &S)
        : radii_0{R}
        , alphas{A}
        , betas{B}
        , sigmas{S}
        , centers{coords} {
    // compute the radii
    for (auto i = 0; i < this->radii_0.size(); ++i) { this->radii.push_back(this->radii_0[i] * this->alphas[i] + this->betas[i] * this->sigmas[i]); }

    auto p_gradcavity = [&cs = this->centers, &rs = this->radii, &ws = this->sigmas](const mrcpp::Coord<3> &r, int index) { return detail::gradCavity(r, index, cs, rs, ws); };
    for (auto i = 0; i < 3; i++) {
        this->gradvector.push_back([i, p_gradcavity](const mrcpp::Coord<3> &r) -> double { return p_gradcavity(r, i); });
    }
}

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
Cavity::Cavity(const std::vector<mrcpp::Coord<3>> &coords, const std::vector<double> &R, double sigma)
        : Cavity(coords, R, std::vector<double>(R.size(), 1.0), std::vector<double>(R.size(), 0.0), std::vector<double>(R.size(), sigma)) {}

/** @brief Evaluates the value of the cavity at a 3D point \f$\mathbf{r}\f$
 *  @param r coordinate of 3D point at which the Cavity is to be evaluated at.
 *  @return double value of the Cavity at point \f$\mathbf{r}\f$
 */
double Cavity::evalf(const mrcpp::Coord<3> &r) const {
    auto C = 1.0;
    for (int i = 0; i < centers.size(); i++) {
        auto center = this->centers[i];
        auto radius = this->radii[i];
        auto sigma = this->sigmas[i];

        auto s = math_utils::calc_distance(center, r) - radius;
        auto Theta = 0.5 * (1 + std::erf(s / sigma));
        auto Ci = 1 - Theta;
        C *= 1 - Ci;
    }
    C = 1 - C;
    return C;
}

bool Cavity::isVisibleAtScale(int scale, int nQuadPts) const {

    auto max_sigma = *std::max_element(this->sigmas.cbegin(), this->sigmas.cend());
    auto visibleScale = static_cast<int>(-std::floor(std::log2(nQuadPts * 2.0 * max_sigma)));

    return (scale < visibleScale);
}

bool Cavity::isZeroOnInterval(const double *a, const double *b) const {
    for (int k = 0; k < this->centers.size(); ++k) {
        auto center = this->centers[k];
        auto radius = this->radii[k];
        auto sigma = this->sigmas[k];
        for (int i = 0; i < 3; ++i) {
            auto cavityMinOut = (center[i] - radius) - 3.0 * sigma;
            auto cavityMinIn = (center[i] - radius) + 3.0 * sigma;
            auto cavityMaxIn = (center[i] + radius) - 3.0 * sigma;
            auto cavityMaxOut = (center[i] + radius) + 3.0 * sigma;

            if (a[i] > cavityMaxOut || (a[i] > cavityMinIn && b[i] < cavityMaxIn) || b[i] < cavityMinOut) { return true; }
        }
    }
    return false;
}
} // namespace mrchem
