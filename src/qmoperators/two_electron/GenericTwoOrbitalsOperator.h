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

#include "GenericTwoOrbitalsPotential.h"
#include "tensor/RankZeroOperator.h"
#include "qmfunctions/Orbital.h"

namespace mrchem {

class GenericTwoOrbitalsOperator: public RankZeroOperator{
public:
    explicit GenericTwoOrbitalsOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi = nullptr, bool mpi_share = false){
        potential = std::make_shared<GenericTwoOrbitalsPotential>(P, Phi, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroOperator &g = (*this);
        g = potential;
        g.name() = "g";
    }
    virtual ~GenericTwoOrbitalsOperator() = default;

    void setup(std::shared_ptr<OrbitalVector> Phi, double prec){
        potential->setup(Phi,prec);
    }

    void set_pair(int j, int l){
        potential->set_pair(j,l);
    }

    Orbital apply(Orbital inp){
        return potential->apply(inp);
    }

private:
    std::shared_ptr<GenericTwoOrbitalsPotential> potential{nullptr};
};

} // namespace mrchem