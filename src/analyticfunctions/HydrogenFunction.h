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

#include "MRCPP/MWFunctions"
#include <array>

namespace mrchem {

class RadialFunction final : public mrcpp::RepresentableFunction<1> {
public:
    RadialFunction(int n, int l, double Z);

    double evalf(const mrcpp::Coord<1> &r) const override;

    friend class HydrogenFunction;

protected:
    const int N;
    const int L;
    double c_0;
    double c_1;

    double calcConstant(double Z) const;
    double evalfPoly(double r) const;

    double calcStdDev() const { return std::pow(2.0 * this->c_1, -0.5); }
    bool isVisibleAtScale(int scale, int nQuadPts) const override;
};

class AngularFunction final : public mrcpp::RepresentableFunction<3> {
public:
    AngularFunction(int l, int m);

    double evalf(const mrcpp::Coord<3> &r) const override;

protected:
    const int L;
    const int M;
    double c_0;

    double calcConstant() const;
    double evalfPoly(const mrcpp::Coord<3> &q) const;
};

class HydrogenFunction final : public mrcpp::RepresentableFunction<3> {
public:
    HydrogenFunction(int n, int l, int m, double Z = 1.0);
    HydrogenFunction(int n, int l, int m, double Z, const mrcpp::Coord<3> &o);

    double evalf(const mrcpp::Coord<3> &r) const override;

protected:
    mrcpp::Coord<3> origin;
    RadialFunction R;
    AngularFunction Y;

    bool isVisibleAtScale(int scale, int nQuadPts) const override { return this->R.isVisibleAtScale(scale, nQuadPts); }
    bool isZeroOnInterval(const double *a, const double *b) const override;
};

} // namespace mrchem
