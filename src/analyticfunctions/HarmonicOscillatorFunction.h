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

#pragma once

#include <array>

#include "MRCPP/MWFunctions"
#include "chemistry/PhysicalConstants.h"

#include "mrchem.h"

namespace mrchem {

class HarmonicOscillator1D final {
public:
    HarmonicOscillator1D(int n, double m, double k, double o)
            : nu(n)
            , alpha(std::pow(m * k, 1.0 / 4.0))
            , origin(o) {}

    double operator()(double x) const {
        double ax = this->alpha * (x - this->origin);
        double ap = std::sqrt(this->alpha / std::sqrt(mrcpp::pi));
        double N_nu = std::sqrt(N2(this->nu));
        return ap * N_nu * H(this->nu, ax) * exp(-ax * ax / 2.0);
    }

protected:
    const double nu;
    const double alpha;
    const double origin;

    double N2(int nu) const {
        if (nu == 0) return 1.0;
        return 1.0 / (2.0 * nu) * N2(nu - 1);
    }
    double H(int nu, double x) const {
        if (nu == 0) return 1.0;
        if (nu == 1) return 2.0 * x;
        return 2.0 * x * H(nu - 1, x) - 2.0 * (nu - 1.0) * H(nu - 2, x);
    }
};

class HarmonicOscillatorFunction final : public mrcpp::RepresentableFunction<3> {
public:
    HarmonicOscillatorFunction(int n[3], double m = 1.0, double *k = nullptr, const mrcpp::Coord<3> &o = {0.0, 0.0, 0.0})
            : fx(n[0], m, ((k != nullptr) ? k[0] : 1.0), o[0])
            , fy(n[1], m, ((k != nullptr) ? k[1] : 1.0), o[1])
            , fz(n[2], m, ((k != nullptr) ? k[2] : 1.0), o[2]) {}

    double evalf(const mrcpp::Coord<3> &r) const override { return fx(r[0]) * fy(r[1]) * fz(r[2]); }

protected:
    const HarmonicOscillator1D fx;
    const HarmonicOscillator1D fy;
    const HarmonicOscillator1D fz;
};

} // namespace mrchem
