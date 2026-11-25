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

namespace mrdft {

/**
 * @class XCLib
 * @brief Abstract exchange–correlation library interface
 * @details This class defines the common interface between MRDFT and the external
 * xc libraries Libxc and XCFun
 */
class XCLib {
public:
    /**
     * @brief Construct a new XC library interface
     * @param[in] spin_enabled If true, configured for
     *                         spin‑polarized calculations, otherwise
     *                         spin‑unpolarized
     */
    explicit XCLib(bool spin_enabled) { this->spin = spin_enabled; }

    virtual ~XCLib() = default;

    XCLib(const XCLib &) = delete;
    XCLib &operator=(const XCLib &) = delete;

    bool spin;                                                    ///< @brief True if density is spin-polarized
    std::vector<std::pair<std::string, double>> functional_specs; ///< @brief Configured functionals and scaling coefficients
    virtual std::string name() const = 0;                         ///< @brief Name of the library (e.g. "XCFun", "Libxc")
    virtual double getAmountExx() const = 0;                      ///< @brief Get total fraction of exact exchange
    virtual int getnOut() const = 0;                              ///< @brief Determine the number of output channels (eg. LDA: 2, spin-polarized LDA: 3, GGA: 5, ...)

    /**
     * @brief Configure a functional object in the backend
     * @param[in] name   XC functional name (e.g. "PBE", "B3LYP",
     *                   Libxc shorthand, or XCFun name)
     * @param[in] c      Scaling coefficient associated with @p name
     */
    virtual void setFunctional(const std::string &name, double c) = 0;

    /**
     * @brief Add a functional specification to #functional_specs
     * @param[in] name XC functional name (e.g. "PBE", "B3LYP",
     *                 Libxc shorthand, or XCFun name)
     * @param[in] c    Scaling coefficient associated with @p name
     * @details
     * This utility is used inside `setFunctional()` implementations
     * to record what has been configured for later reporting
     */
    void addFunctionalSpec(const std::string &name, double c) { functional_specs.emplace_back(name, c); }

    /**
     * @brief Inspect functional families and sets evaluation variables for the library
     * @param[out] lda  True if any functional is LDA
     * @param[out] gga  True if any functional is GGA
     * @param[out] mgga True if any functional is metaGGA
     * @param[in]  order Order of functional derivatives requested (Libxc supports only order <= 1)
     * @param[in]  gamma If true, gamma‑type derivatives are requested (unused for Libxc)
     */
    virtual void initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, int order, bool gamma) = 0;

    /**
     * @brief Evaluate configured functionals on a batch of grid points
     * @param[in]  inp    Input matrix with one column per grid point and
     *                    one row per input component
     * @param[out] out    Output matrix with energy densities and potentials
     * @param[in]  nPts   Number of grid points (cols in @p inp/@p out)
     */
    virtual void callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts) const = 0;

    /**
     * @brief Print information and references for all active xc functionals
     * @param[in] out_txt_width Max width for output used to break lines
     */
    virtual void printFunctionalReference(int out_txt_width) const = 0;

    /**
     * @brief Set the density cutoff
     * @param[in] cutoff Density threshold. XCFun skips grid points with density below the cutoff.
     *                   Libxc uses it as its density threshold (::xc_func_set_dens_threshold);
     *                   values <= 0 mean Libxc's built-in default
     */
    virtual void setCutoff(double cutoff) = 0;

    /**
     * @brief Print a common xc functional reference header
     */
    static void printReferenceHeader();

    /**
     * @brief Helper function for word‑wrapping reference strings
     * @param[in] str       String to be wrapped
     * @param[in] txt_width Maximum line width
     * @param[in] indent    Number of spaces to indent each wrapped line, used for tabbed lines
     */
    static void printWrap(const std::string &str, std::size_t txt_width, int indent = 0);

protected:
    double cutoff{-1.0};
};

} // namespace mrdft