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

#include "tensor/RankZeroOperator.h"

namespace mrchem {

class NuclearFunction;

class NuclearOperator final : public RankZeroOperator {
public:
    NuclearOperator(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false, const std::string &model = "point_like");

    inline void add(const NuclearOperator &op) {
        this->V_func.add(1.0, op.V_func);
    }

private:
    mrcpp::ComplexFunction V_func; // The MW Function representation of the potential

    void setupLocalPotential(NuclearFunction &f_loc, const Nuclei &nucs, double smooth_prec) const;
    void allreducePotential(double prec, mrcpp::ComplexFunction &V_tot, mrcpp::ComplexFunction &V_loc) const;
};

} // namespace mrchem
