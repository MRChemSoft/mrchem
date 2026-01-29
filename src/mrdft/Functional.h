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
#include <iostream>
#include <set>

#include <Eigen/Core>
#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <XCFun/xcfun.h>
#include <xc_funcs.h>
#include <xc.h>

namespace mrdft {

using XC_p = std::unique_ptr<xcfun_t, decltype(&xcfun_delete)>;

/**
 * @class Functional
 * @brief Abstract base class for Exchange-Correlation functionals
 * @details This class provides the interface for evaluating XC potentials
 * on the Multi-Resolution Analysis (MRA) grid using either Libxc or XCFun
 */
class Functional {
public:
    /**
     * @brief Constructor for the Functional base class
     * @param[in] k The polynomial order for the MRA basis
     * @param[in] f The XCFun handle (ownership is transferred)
     */
    Functional(int k, XC_p &f)
            : order(k)
            , xcfun(std::move(f)) {}
    virtual ~Functional() = default;

    /** @brief  Evaluates XC functional and derivatives for a given NodeIndex
     *
     * The electronic densities (total/alpha/beta) are given as input.
     * The values of the zero order densities and their gradient are sent to xcfun.
     * The output of xcfun must then be combined ("contract") with the gradients
     * of the higher order densities.
     *
     * XCFunctional output (with k=1 and explicit derivatives):
     *
     * LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho}\right) \f$
     *
     * GGA: \f$ \left(F_{xc},
     *  \frac{\partial F_{xc}}{\partial \rho},
     *  \frac{\partial F_{xc}}{\partial \rho_x},
     *  \frac{\partial F_{xc}}{\partial \rho_y},
     *  \frac{\partial F_{xc}}{\partial \rho_z}\right) \f$
     *
     * Spin LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho^\alpha},
     *  \frac{\partial F_{xc}}{\partial \rho^\beta}\right) \f$
     *
     * Spin GGA: \f$ \left(F_{xc},
     *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
     *  \frac{\partial F_{xc}}{\partial \rho^\beta},
     *  \frac{\partial F_{xc}}{\partial \rho_x^\alpha},
     *  \frac{\partial F_{xc}}{\partial \rho_y^\alpha},
     *  \frac{\partial F_{xc}}{\partial \rho_z^\alpha},
     *  \frac{\partial F_{xc}}{\partial \rho_x^\beta},
     *  \frac{\partial F_{xc}}{\partial \rho_y^\beta},
     *  \frac{\partial F_{xc}}{\partial \rho_z^\beta}
     *  \right) \f$
     *
     * XCFunctional output (with k=1 and gamma-type derivatives):
     *
     * GGA: \f$ \left(F_{xc},
     *  \frac{\partial F_{xc}}{\partial \rho},
     *  \frac{\partial F_{xc}}{\partial \gamma} \f$
     *
     * Spin GGA: \f$ \left(F_{xc},
     *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
     *  \frac{\partial F_{xc}}{\partial \rho^\beta },
     *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \alpha}},
     *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \beta }},
     *  \frac{\partial F_{xc}}{\partial \gamma^{\beta  \beta }}
     *  \right) \f$
     *
     * @param[in] inp Input values
     * @param[out] xcNodes Output values
     *
     */
    void makepot(mrcpp::FunctionTreeVector<3> &inp, std::vector<mrcpp::FunctionNode<3> *> xcNodes) const;

    /**
     * Setters
     */
    void setLogGradient(bool log) { log_grad = log; }    ///< @brief Set whether to use logarithmic gradient transformations
    void setDensityCutoff(double cut) { cutoff = cut; }  ///< @brief Set the density threshold below which density is set to 0
    void setDerivOp(std::unique_ptr<mrcpp::DerivativeOperator<3>> &d) { derivOp = std::move(d); }   ///< @brief Set the numerical derivative operator for gradient-based functionals

    /**
     * Functional type querying
     */
    bool isLDA() const { return not (isGGA() or isMetaGGA()); }          ///< @return True if functional is LDA type (not a GGA or meta-GGA)
    bool isHybrid() const { return (std::abs(amountEXX()) > 1.0e-10); }  ///< @return True if functional is a hybrid (includes exact exchange)
    virtual bool isSpin() const = 0;                                     ///< @brief Returns True if the functional object is spin-polarized
    virtual bool isGGA() const = 0;                                      ///< @brief Returns True if the functional is a GGA
    virtual bool isMetaGGA() const = 0;                                  ///< @brief Returns True if the functional is a Meta-GGA

    virtual int numIn() const = 0;      ///< Fetches number of variables in the input matrix
    virtual int numOut() const = 0;     ///< Fetches number of variables in the output matrix

    /**
     * @brief Fetches the amount of exact exchange needed for a given functional
     * @return The total fraction of exx to be added to the functional
     */
    double amountEXX() const;
    double XCenergy = 0.0;          ///< @brief Stores calculated xc energy for the current state

    /**
     * @brief Evaluates the functional on a set of grid points
     * @param[in] inp  Matrix of input values (density, gradient,...)
     * @param[out] out out_data Matrix of output values (energy, potential, ...)
     * @details Each column corresponds to one grid point
     * From a performance point of view, (in pre and postprocessing) it is much more
     * efficient to have the two consecutive points in two consecutive adresses in memory
     */
    Eigen::MatrixXd evaluate(Eigen::MatrixXd &inp) const;

    /**
     * @brief Evaluates the functional on a set of grid points. Transposed version of Functional::evaluate()
     * @param[in] inp Matrix of input values (density, gradient,...)
     * @param[out] out_trans Matrix of output values (energy, potential, ...)
     * @details Each row corresponds to one grid point
     */
    Eigen::MatrixXd evaluate_transposed(Eigen::MatrixXd &inp) const;
    
    bool libxc;                                 ///< @brief Flag indicating if Libxc is active (True if "DFT {xc_library = libxc}" in input file)
    std::vector<xc_func_type> libxc_objects;    ///< @brief Vector of initialized Libxc functionals
    std::vector<double> libxc_coefs;            ///< @brief Vector scaling coefficients for each functional in libxc_objects
    
    /**
     * @brief Prints the splash screens, version info, and references for the 
     * active xc libraries and functionals
     * @details If Libxc is used, it iterates through 
     * all initialized functional objects to print their specific DOIs
     */
    void print_functional_references() const;
    
    /**
     * @brief Transfers ownership of Libxc functional objects and their scaling 
     * coefficients to the Functional instance
     * @param[in] libxc_objects_ Vector of initialized Libxc functionals
     * @param[in] libxc_coefs_   Vector of corresponding weights of the initialized Libxc functionals
     */
    void set_libxc_functional_object(std::vector<xc_func_type> libxc_objects_, std::vector<double> libxc_coefs_);

    friend class MRDFT;

protected:
    const int order;            ///< @brief Order of contraction ????
    bool log_grad{false};       ///< @brief Toggle for logarithmic gradient
    double cutoff{-1.0};        ///< @brief Density threshold
    Eigen::VectorXi d_mask;     ///< @brief 
    Eigen::MatrixXi xc_mask;    ///< @brief 
    XC_p xcfun;                 ///< XCFun library handle
    std::unique_ptr<mrcpp::DerivativeOperator<3>> derivOp{nullptr};  ///< @brief Operator used to compute gradients

    /**
     * @brief Run a collection of grid points through Libxc or XCFun
     * @param[in] inp  Matrix of input values, where each row is one grid point
     * @param[out] out Matrix of output values
     */
    void evaluate_data(const Eigen::MatrixXd & inp, Eigen::MatrixXd &out) const;

    /**
     * @brief Contracts a collection of grid points
     * @details This is used to implement the chain rule for functionals involving 
     * gradients or when calculating higher-order properties
     * @param[in] xc_data xc_data Matrix of functional partial derivative values
     * @param[in] d_data  d_data Matrix of density input values
     * @param[out] out_data Matrix of contracted output values
     * @details Each row corresponds to one grid point
     */
    Eigen::MatrixXd contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;
    /**
     * @brief Contracts a collection of grid points. Transposed version of Functional::contract
     * @details This is used to implement the chain rule for functionals involving 
     * gradients or when calculating higher-order properties
     * @param[in] xc_data xc_data Matrix of functional partial derivative values
     * @param[in] d_data  d_data Matrix of density input values
     * @param[out] out_data Matrix of contracted output values
     * @details Each column corresponds to one grid point
     */
    Eigen::MatrixXd contract_transposed(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;

    virtual int getCtrInputLength() const = 0;                          ///< @brief Expected number of input components for the contraction step
    virtual int getCtrOutputLength() const = 0;                         ///< @brief Expected number of output components for the contraction step
    virtual void clear() = 0;                                           ///< @brief Clears internal functions
    virtual mrcpp::FunctionTreeVector<3> setupXCInput() = 0;            ///< @brief Configures input for evaluation
    virtual mrcpp::FunctionTreeVector<3> setupCtrInput() = 0;           ///< @brief Configures input for contraction
    virtual void preprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;     ///< @brief Collects input functions for evaluation
    virtual mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) = 0; ///< @brief Computes final output functions

};

} // namespace mrdft
