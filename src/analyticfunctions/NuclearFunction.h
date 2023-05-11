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
#include <vector>

#include <MRCPP/MWFunctions>

#include "chemistry/Nucleus.h"

namespace mrchem {

class NuclearFunction final : public mrcpp::RepresentableFunction<3> {
public:
    NuclearFunction(){};
    NuclearFunction(const Nuclei &nucs, double smooth_prec = -1.0, double prec = -1.0);

    double evalf(const mrcpp::Coord<3> &r) const override;

    void push_back(const std::string &atom, const mrcpp::Coord<3> &r, double c);
    void push_back(const Nucleus &nuc, double c);

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

    double getPrec() { return this->prec; }

    bool isVisibleAtScale(int scale, int nQuadPts) const override;
    bool isZeroOnInterval(const double *a, const double *b) const override;

protected:
    Nuclei nuclei;
    double prec;
    std::vector<double> smooth;
    std::vector<double> minPot;
    std::string smooth_method = "very_smooth"; // can be: parabola, minimum or very_smooth
};

namespace detail {
/*! @brief Compute nucleus- and precision-dependent smoothing parameter */
inline auto nuclear_potential_smoothing(double smooth_prec, double Z) -> double {
    double tmp = 0.00435 * smooth_prec / std::pow(Z, 5.0);
    return std::cbrt(tmp);
}
} // namespace detail

} // namespace mrchem
