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

#include "xc.h"
#include <vector>
#include <Eigen/Core>
#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/trees/FunctionTreeVector.h>
#include <MRCPP/trees/FunctionTree.h>
#include "xc_utils.h"

#include "Factory.h"
#include "LibXC.h"

namespace mrdft {

class libxc_LDA final : public LibXC {
public:
    libxc_LDA(int order, std::vector<mrlibxc::LibXCData> &functionals);
    ~libxc_LDA() override = default;


    bool isSpin() const override { return false; }

private:
    mrcpp::FunctionTree<3> rho_ft;
    XC_p xc_;
    std::vector<mrlibxc::LibXCData> functionals_;
    mrcpp::FunctionTreeVector<3> rho;
    int order_;

    int getCtrInputLength() const override;
    int getCtrOutputLength() const override { return 2; }

    void clear() override;
    mrcpp::FunctionTreeVector<3> setupXCInput() override;
    mrcpp::FunctionTreeVector<3> setupCtrInput() override;

    void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
    mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;

};

} // namespace mrdft

