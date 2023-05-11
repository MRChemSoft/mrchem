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

#include <cmath>

#include <MRCPP/MWFunctions>

namespace mrchem {

class NuclearGradientFunction : public mrcpp::RepresentableFunction<3> {
public:
    /*!
     *  @brief NuclearGradientFunction represents the function: Z * [x,y,z]/|r - o|^3
     *  @param d: Cartesian component to evaluate (0: x, 1: y, 2: z)
     *  @param z: Nuclear charge of nucleus
     *  @param o: Coordinate of origin
     *  @param c: Nuclei- and precision-dependent smoothing parameter
     */
    NuclearGradientFunction(int d, double z, const mrcpp::Coord<3> &o, double c)
            : D(d)
            , C(c)
            , Z(z)
            , R(o) {}

    double evalf(const mrcpp::Coord<3> &r) const override;

    bool isVisibleAtScale(int scale, int nQuadPts) const override;
    bool isZeroOnInterval(const double *a, const double *b) const override;

protected:
    int D;             ///< Cartesian direction
    double C;          ///< Smoothing parameter (nuclei- and precision-dependent)
    double Z;          ///< Nuclear charge
    mrcpp::Coord<3> R; ///< Nuclear coordinate

    double du_dr(double r1) const;
};

namespace detail {
/*! @brief Compute nucleus- and precision-dependent smoothing parameter */
inline auto nuclear_gradient_smoothing(double prec, double Z, int N) -> double {
    auto tmp = prec / (0.00435 * std::pow(Z, 5) * N);
    return std::cbrt(tmp);
}
} // namespace detail

} // namespace mrchem
