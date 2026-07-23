// /*
//  * MRChem, a numerical real-space code for molecular electronic structure
//  * calculations within the self-consistent field (SCF) approximations of quantum
//  * chemistry (Hartree-Fock and Density Functional Theory).
//  * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
//  *
//  * This file is part of MRChem.
//  *
//  * MRChem is free software: you can redistribute it and/or modify
//  * it under the terms of the GNU Lesser General Public License as published by
//  * the Free Software Foundation, either version 3 of the License, or
//  * (at your option) any later version.
//  *
//  * MRChem is distributed in the hope that it will be useful,
//  * but WITHOUT ANY WARRANTY; without even the implied warranty of
//  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  * GNU Lesser General Public License for more details.
//  *
//  * You should have received a copy of the GNU Lesser General Public License
//  * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
//  *
//  * For information on the complete list of contributors to MRChem, see:
//  * <https://mrchem.readthedocs.io/>
//  */

// #pragma once

// #include <XCFun/xcfun.h>

// #include "Functional.h"
// #include "Factory.h" // only to call Factory::libxc

// namespace mrdft {

// class SpinmGGA final : public Functional {
// public:
//     SpinmGGA(int k, XC_p &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> &d);
//     ~SpinmGGA() override = default;

//     bool isSpin() const override { return true; }
//     bool isGGA() const override { return true; }
//     bool isMetaGGA() const override { return true; }
//     int numIn() const override { return 10; }  // rho_a, rho_b, grad_a×3, grad_b×3, tau_a, tau_b
//     int numOut() const override { 
//         if (Factory::libxc) {
//             // f_xc, v_rho_a, v_rho_b, grad_a×3, grad_b×3, v_tau_a, v_tau_b
//             return 11;
//         } else {
//             return xcfun_output_length(xcfun.get());
//         }
//     }

// private:
//     std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
//     mrcpp::FunctionTreeVector<3> rho_a;
//     mrcpp::FunctionTreeVector<3> rho_b;
//     mrcpp::FunctionTreeVector<3> grad_a;
//     mrcpp::FunctionTreeVector<3> grad_b;
//     mrcpp::FunctionTreeVector<3> tau_a;
//     mrcpp::FunctionTreeVector<3> tau_b;

//     int getCtrInputLength() const override;
//     int getCtrOutputLength() const override { return 11; }

//     void clear() override;
//     virtual mrcpp::FunctionTreeVector<3> setupXCInput() override;
//     virtual mrcpp::FunctionTreeVector<3> setupCtrInput() override;

//     void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
//     mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;
// };

// } // namespace mrdft
