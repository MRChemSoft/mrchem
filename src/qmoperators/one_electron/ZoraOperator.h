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
#include <string>

namespace mrchem {

class QMPotential;

/**
 * @class ZoraOperator
 * @brief Implements chi = kappa - 1 relativistic dampening function. This has to be done in order to
 * avoid numerical instabilities in the ZORA operator. Whenever this operator is applied, the + 1 has to be added to chi i
 * in order to get the kappa operator. kappa * phi = chi * phi + phi
 * This has to be done manually.
 */
class ZoraOperator final : public RankZeroOperator {
public:
    /**
     * @brief Constructor for the ZoraOperator that contains the chi = kappa - 1 function.
     * @param vz The potential used to calculate the kappa function.
     * @param c Speed of light.
     * @param proj_prec The precision of the MW projection.
     * @param inverse If true, the inverse of the chi function is calculated.
     */
    ZoraOperator(QMPotential &vz, double c, double proj_prec, bool inverse = false);

    ZoraOperator(std::shared_ptr<QMPotential> relativisticDampening, std::string name);

};

} // namespace mrchem
