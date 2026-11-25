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

#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <memory>
#include "xc_libraries/XCLib.h"

namespace mrdft {

using XCLib_p = std::unique_ptr<XCLib>;

/**
 * @class Functional
 * @brief Abstract base class for Exchange-Correlation functionals
 * @details Provides the interface for evaluating XC potentials
 * on the Multi-Resolution Analysis (MRA) grid using either Libxc or XCFun
 */
class Functional {
public:
    /**
     * @brief Constructor for the Functional base class
     * @param[in] k Order of the functional derivatives (0: energy, 1: potential, etc.)
     * @param[in] s True for a spin-polarized functional (alpha and beta densities), false for spin-restricted
     * @param[in] f The XC library backend (ownership is transferred)
     */
    Functional(int k, bool s, XCLib_p &f)
            : order(k)
            , spin{s} // braces: Orbital.h defines a function-like macro spin()
            , xclib(std::move(f)) {}
    virtual ~Functional() = default;

    /** @brief  Evaluates XC functional and derivatives for a given NodeIndex
     *
     * The electronic densities (total/alpha/beta) are given as input.
     * The values of the zero order densities and their gradient are sent to the XC library.
     * The output of the XC library must then be combined ("contract") with the gradients
     * of the higher order densities.
     *
     * XC library output (with k=1 and explicit derivatives):
     *
     * LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho}\right) \f$
     *
     * GGA: \f$ \left(F_{xc},
     *  \frac{\partial F_{xc}}{\partial \rho},
     *  \frac{\partial F_{xc}}{\partial \rho_x},
     *  \frac{\partial F_{xc}}{\partial \rho_y},
     *  \frac{\partial F_{xc}}{\partial \rho_z}\right) \f$
     *
     * LDA, spin-polarized: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho^\alpha},
     *  \frac{\partial F_{xc}}{\partial \rho^\beta}\right) \f$
     *
     * GGA, spin-polarized: \f$ \left(F_{xc},
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
     * XC library output (with k=1 and gamma-type derivatives):
     *
     * GGA: \f$ \left(F_{xc},
     *  \frac{\partial F_{xc}}{\partial \rho},
     *  \frac{\partial F_{xc}}{\partial \gamma} \f$
     *
     * GGA, spin-polarized: \f$ \left(F_{xc},
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
    void setLogGradient(bool log) { log_grad = log; }                                             ///< @brief Set whether to use logarithmic gradient transformations (currently unused, see log_grad)
    void setDensityCutoff(double cut) { cutoff = cut; }                                           ///< @brief Set the density threshold below which density is set to 0
    void setDerivOp(const std::shared_ptr<mrcpp::DerivativeOperator<3>> &d) { derivOp = d; }        ///< @brief Set the numerical derivative operator for gradient-based functionals

    /**
     * Functional type querying
     */
    bool isLDA() const { return not(isGGA() or isMetaGGA()); }          ///< @return True if functional is LDA type (not a GGA or meta-GGA)
    bool isHybrid() const { return (std::abs(amountEXX()) > 1.0e-10); } ///< @return True if functional is a hybrid (includes exact exchange)
    bool isSpin() const { return spin; }                                ///< @brief Returns True if the functional object is spin-polarized
    virtual bool isGGA() const = 0;                                     ///< @brief Returns True if the functional is a GGA
    virtual bool isMetaGGA() const = 0;                                 ///< @brief Returns True if the functional is a Meta-GGA

    int numIn() const { return densityChannels() * (usesGradients() ? 4 : 1); } ///< Fetches number of variables in the input matrix (per density channel: 1 density, +3 gradient components if usesGradients())
    virtual int numOut() const = 0; ///< Fetches number of variables in the output matrix

    /**
     * @brief Fetches the amount of exact exchange needed for a given functional
     * @return The total fraction of exx to be added to the functional
     */
    double amountEXX() const;
    double XCenergy = 0.0; ///< @brief Stores calculated xc energy for the current state

    /**
     * @brief Evaluates the functional on a set of grid points
     * @param[in] inp Matrix of input values (density, gradient,...)
     * @return Matrix of output values (energy, potential, ...)
     * @details Each row corresponds to one grid point Previously called evaluate_transposed()
     */
    Eigen::MatrixXd evaluate(Eigen::MatrixXd &inp) const;

    /**
     * @brief Prints the splash screens, version info, and references for the
     * active xc libraries and functionals
     * @details If Libxc is used, it iterates through
     * all initialized functional objects to print their specific DOIs
     */
    void print_functional_references() const;

    friend class MRDFT;

protected:
    const int order;                                                ///< @brief Order of functional derivatives. Eg. 0 (energy), 1 (potential), 2 (pot. gradient), etc.
    const bool spin;                                                ///< @brief True if the functional is spin-polarized (alpha and beta densities)
    /** @brief Toggle for logarithmic gradients (xc_utils::log_gradient)
     *
     * Currently unused: no input keyword sets it, so it is always false. makepot() computes the density
     * gradients node by node, while xc_utils::log_gradient works on whole FunctionTrees (log of the full
     * density, then a tree-wide gradient), so it cannot be called from makepot(). If set, makepot() reports
     * an error. Kept so a tree-level gradient path can use it later.
     */
    bool log_grad{false};
    double cutoff{-1.0};                                            ///< @brief Density threshold. Currently unused (the cutoff is applied in XCLib); used by the meta-GGA work on the mggas branch
    Eigen::VectorXi d_mask;                                         ///< @brief density and derivative(s) mask vector for response calculations
    Eigen::MatrixXi xc_mask;                                        ///< @brief functional and derivative(s) mask vector for response calculations
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivOp{nullptr}; ///< @brief Operator used to compute gradients
    std::unique_ptr<XCLib> xclib;                                   ///< @brief Handle for exchange–correlation library interface

    int densityChannels() const { return spin ? 2 : 1; } ///< @brief Returns number of density channels (spin-restricted -> 1, spin-polarized -> 2)
    virtual bool usesGradients() const = 0;  ///< @brief Returns true if functional type uses gradients (eg. true for GGAs or metaGGAs)

    /**
     * @brief Run a collection of grid points through Libxc or XCFun
     * @param[in] inp  Matrix of input values, where each row is one grid point
     * @param[out] out Matrix of output values
     */
    void evaluate_data(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;

    /**
     * @brief Contracts a collection of grid points
     * @details This is used to implement the chain rule for functionals involving
     * gradients or when calculating higher-order properties
     * @param[in] xc_data xc_data Matrix of functional partial derivative values
     * @param[in] d_data  d_data Matrix of density input values
     * @return Matrix of contracted output values. Each row corresponds to one grid point
     */
    Eigen::MatrixXd contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;

    int getCtrOutputLength() const { return numIn() + 1; } ///< @brief Expected number of output components for the contraction step (energy + one per input component)
};

} // namespace mrdft
