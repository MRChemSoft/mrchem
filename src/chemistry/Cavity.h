/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/MWFunctions"
#include <array>
#include <vector>

using namespace mrcpp;

namespace mrchem {


class Cavity : public mrcpp::RepresentableFunction<3> {
public:
  Cavity(std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope);
  Cavity(std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope, double e_i, double e_o);
  void eval_epsilon(bool argument, bool implement);
  double evalf(const double *r) const override;
  double evalf(const std::array<double, 3> &r) const {return evalf(r.data());};

protected:
    double e_i = 1;
    double e_o = 2;
    bool is_eps = false;
    bool is_linear = false;
    std::vector<mrcpp::Coord<3>> pos;
    std::vector<double> R;
    double d;
    double dcoeff;


} // namespace mrchem
