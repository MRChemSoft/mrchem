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
#include <XCFun/xcfun.h>
#include <xc_funcs.h>
#include <xc.h>

#include "MRDFT.h"

namespace mrdft {

/**
 * @class Factory
 * @brief Building class for MRDFT objects
 * @details Manages the different xc libraries, mapping of functional names and builds the functional objects 
 * with the required parameters to initialize a DFT calculation
 */
class Factory final {
public:
    /**
     * @brief Construct a new Factory object. Initializes the MRA reference and creates the XCFun handle
     * @param[in] MRA The Multi-Resolution Analysis object providing the grid 
     * and basis functions for the calculation
     */
    Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA);

    ~Factory() = default;   ///< @brief Default destructor

    /*
     * Setters
     */
    void setSpin(bool s) { spin = s; }                        ///< Set spin polarization (true for unrestricted/spin-polarized) */
    void setOrder(int k) { order = k; }                       ///< Set the polynomial order for the MRA basis
    void setUseGamma(bool g) { gamma = g; }                   ///< Toggle between gamma-type and explicit derivatives
    void setLogGradient(bool lg) { log_grad = lg; }           ///< Toggle the use of logarithmic gradients
    void setDensityCutoff(double c) { cutoff = c; }           ///< Set the threshold for neglecting low-density regions
    void setLibxc(bool libxc_) {libxc = libxc_; }             ///< Toggle between Libxc (true) and XCFun (false) backends
    void setDerivative(const std::string &n) { diff_s = n; }  ///< Set derivative operator type (e.g., "bspline", "abgv_00")

    /**
     * @brief Configures the xc functional
     * 
     * @param[in] name The name of the xc functional (e.g., "PBE", "B3LYP")
     * @param[in] c    A global scaling coefficient applied to the functional.
     * @throws MSG_ABORT If a mapped Libxc ID is incompatible with the linked Libxc version
     * @note Depending on the chosen library, this method either initializes
     * one or more Libxc functional objects using the MapFuncName function or sets
     * the functional parameters in the XCFun backend
     */
    void setFunctional(const std::string &n, double c = 1.0);

    /**
     * @brief Build a MRDFT object from the currently defined parameters
     * @details Performs the following steps:
     * 1.  **Grid Initialization**: Creates a multi-resolution grid based on the MRA
     * 2.  **Library dependent initiation**: If using Libxc, iterates through functional objects
     * to ensure they belong to supported families (LDA/GGA, not meta-GGA or range separated). 
     * If using XCFun, sets evaluation parameters, mode and order
     * 3.  **Operator Selection**: Assigns numerical derivative operators (BSpline or ABGV)
     * required for GGAs
     * 4.  **Functional Instantiation**: Selects the appropriate concrete implementation
     * (SpinLDA, SpinGGA, LDA, or GGA) based on spin and gradient requirements.
     * 5.  **State Sync**: Passes functional objects, density cutoffs, and derivative 
     * schemes to the functional
     * @return std::unique_ptr<MRDFT> A pointer to the assembled Multi-Resolution DFT object.
     * @throws MSG_ABORT If unsupported functional families are detected in the Libxc case
     * (eg. meta-GGAs and range separated functionals)
     */
    std::unique_ptr<MRDFT> build();

    static bool libxc;      ///< @brief Flag indicating if Libxc is active (True if "DFT {xc_library = libxc}" in input file)

private:
private:
    int order{1};                  ///< Polynomial order of the Multi-Resolution Analysis (MRA) basis
    bool spin{false};              ///< If true, perform unrestricted calculations
    bool gamma{false};             ///< If true, use gamma-type derivatives (gradient squared) instead of explicit components
    bool log_grad{false};          ///< Toggle for using logarithmic gradient transformations
    double cutoff{-1.0};           ///< Density threshold; values below this are sat to 0
    std::string diff_s{"abgv_00"}; ///< String identifier for the derivative operator type (e.g., "bspline", "abgv_55")

    XC_p xcfun_p;                                         ///< @brief Pointer to the XCFun library handle
    std::unique_ptr<mrcpp::DerivativeOperator<3>> diff_p; ///< @brief Pointer to the numerical derivative operator used for GGA gradients


    std::vector<xc_func_type> libxc_objects;        ///< @brief Vector of initialized Libxc functionals
    std::vector<double> libxc_coefs;                ///< @brief Vector scaling coefficients for each functional in libxc_objects

    /**
     * @brief Maps a functional name string (e.g., "PBE0", "LDA" or "XC_LDA_X", XC_GGA_X_B88) 
     * to its corresponding Libxc IDs and scaling coefficients
     * @note The input `name` is transformed to uppercase internally, making the
     * search case-insensitive
     * @param[in] name    Name of the functional
     * @param[in] ids     Vector to be populated with the IDs used by Libxc
     * @param[in] coefs   Vector to be populated with the corresponding scaling coefficients
     * @throws MSG_ABORT If the name is not a recognized internal shorthand and
     * is not found within the Libxc library
     * @example
     * std::vector<int> ids;
     * std::vector<double> coefs;
     * MapFuncName("LDA", ids, coefs); 
     * // ids: {XC_LDA_C_VWN, XC_LDA_X}, coefs: {1.0, 1.0}
     */
    std::vector<int> mapFunctionalName(const std::string &name) const;
};

} // namespace mrdft
