/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/MWOperators"

#include "tensor/RankZeroOperator.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "qmoperators/one_electron/NablaOperator.h"

namespace mrchem {

class ZoraOperator final : public RankZeroOperator {
public:
    ZoraOperator(QMPotential &V, double c, std::shared_ptr<mrcpp::DerivativeOperator<3>> D, bool mpi_share = false) 
        : light_speed(c)
        , derivative(D) {
        auto V_zora = std::make_shared<QMPotential>(1, mpi_share);
        qmfunction::deep_copy(*V_zora, V);
        
        // Re-map the ZORA function on the same grid as the input potential
        double zfac = 2.0 * c * c;
        auto zmap = [zfac](double val) -> double { return zfac / (zfac - val); };
        V_zora->real().map(zmap);
        this->potential = V_zora;

        // Invoke operator= to assign *this operator
        RankZeroOperator &O = (*this);
        O = V_zora;
        O.name() = "V_zora";
};

public:
    double light_speed;
    std::shared_ptr<QMPotential> potential = std::make_shared<QMPotential>(1, false);
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;

public:
    std::shared_ptr<RankZeroOperator> divKappaOverSqKappa() {
        RankZeroOperator A(this->invSqKappa());
        RankZeroOperator B(this->divSqKappa());
        RankZeroOperator prod = A * B;
        auto out = std::make_shared<RankZeroOperator>(prod);
        return out;
        }

    std::shared_ptr<QMPotential> sqKappa() {
        auto srmap = [](double val) -> double {return std::sqrt(val); };
        auto sq_V_zora = std::make_shared<QMPotential>(1, false);
        qmfunction::deep_copy(*sq_V_zora, *(this->potential));
        sq_V_zora->real().map(srmap);
        return sq_V_zora;
    }

    std::shared_ptr<QMPotential> invSqKappa() {
        auto srmap = [](double val) -> double {return 1 / std::sqrt(val); };
        auto inv_sq_V_zora = std::make_shared<QMPotential>(1, false);
        qmfunction::deep_copy(*inv_sq_V_zora, *(this->potential));
        inv_sq_V_zora->real().map(srmap);
        return inv_sq_V_zora;
    }
 
    std::shared_ptr<QMPotential> divSqKappa() {
        mrcpp::FunctionTreeVector<3> gradK = mrcpp::gradient(*(this->derivative), this->sqKappa()->real());
        auto divK = std::make_shared<QMPotential>(1, false);
        divK->alloc(NUMBER::Real);
        mrcpp::divergence(divK->real(), *(this->derivative), gradK);
        return divK;
    }
};
} // namespace mrchem
