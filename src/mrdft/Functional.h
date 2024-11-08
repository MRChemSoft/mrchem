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

#include <memory>

#include <Eigen/Core>
#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <XCFun/xcfun.h>

namespace mrdft {

using XC_p = std::unique_ptr<xcfun_t, decltype(&xcfun_delete)>;

class Functional {
public:
    Functional(int k, XC_p &f)
            : order(k)
            , xcfun(std::move(f)) {}
    virtual ~Functional() = default;

    void makepot(mrcpp::FunctionTreeVector<3> &inp, std::vector<mrcpp::FunctionNode<3> *> xcNodes) const;

    void setLogGradient(bool log) { log_grad = log; }
    void setDensityCutoff(double cut) { cutoff = cut; }
    void setDerivOp(std::unique_ptr<mrcpp::DerivativeOperator<3>> &d) {derivOp = std::move(d);}

    virtual bool isSpin() const = 0;
    bool isLDA() const { return (not(isGGA() or isMetaGGA())); }
    bool isGGA() const { return xcfun_is_gga(xcfun.get()); }
    bool isMetaGGA() const { return xcfun_is_metagga(xcfun.get()); }
    bool isHybrid() const { return (std::abs(amountEXX()) > 1.0e-10); }
    double amountEXX() const {
        double exx = 0.0;
        xcfun_get(xcfun.get(), "exx", &exx);
        return exx;
    }
    double XCenergy = 0.0;

    Eigen::MatrixXd evaluate(Eigen::MatrixXd &inp) const;
    Eigen::MatrixXd evaluate_transposed(Eigen::MatrixXd &inp) const;
    friend class MRDFT;

protected:
    const int order;
    bool log_grad{false};
    double cutoff{-1.0};
    Eigen::VectorXi d_mask;
    Eigen::MatrixXi xc_mask;
    XC_p xcfun;
    std::unique_ptr<mrcpp::DerivativeOperator<3>> derivOp{nullptr};

    int getXCInputLength() const { return xcfun_input_length(xcfun.get()); }
    int getXCOutputLength() const { return xcfun_output_length(xcfun.get()); }
    virtual int getCtrInputLength() const = 0;
    virtual int getCtrOutputLength() const = 0;

    Eigen::MatrixXd contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;
    Eigen::MatrixXd contract_transposed(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;

    virtual void clear() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupXCInput() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupCtrInput() = 0;

    virtual void preprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;
    virtual mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;
};

} // namespace mrdft
