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

#include "XCLib.h"
#include <xc.h>
#include <xc_funcs.h>

namespace mrdft {

/**
 * @class Libxc
 * @brief Class for Exchange-Correlation functionals using Libxc
 * @details Provides the interface for evaluating xc potentials
 * on the Multi-Resolution Analysis (MRA) grid using Libxc
 */
class Libxc final : public XCLib {

public:
    explicit Libxc(bool spin_enabled)
            : XCLib(spin_enabled) {}

    ~Libxc() override;

    std::vector<xc_func_type *> libxc_objects; ///< @brief Vector of initialized Libxc functionals
    std::vector<double> libxc_coefs;           ///< @brief Vector scaling coefficients for each functional in libxc_objects
    double getAmountExx() const override;
    double getCustomExx() const { return customExx; }
    int getnOut() override;

    void setFunctional(const std::string &name, double c) override;
    void initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, int order, bool gamma) override;
    void callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts) const override;
    void printFunctionalReference(int out_txt_width) const override;

    void setCutoff(double cutoff) override;

    double customExx{0.0};                             ///< @brief Additional exact exchange provided by mapFunctionalName for custom functionals
    void setCustomExx(double exx) { customExx = exx; } ///< @brief Used in mapFunctionalName() to set a custom fraction of exact exchange
    /**
     * @brief Maps a functional name string (eg., "PBE0", "LDA" or "XC_LDA_X", XC_GGA_X_B88) to its corresponding Libxc IDs and scaling coefficients
     * @param[in] name        Name of the functional
     * @param[in] ids         Vector to be populated with the IDs used by Libxc
     * @param[in] coefs       Vector to be populated with the corresponding scaling coefficients
     * @param[in] customExx   Variable for a customized exact exchange sat by a pre-defined functional
     */
    void mapFunctionalName(std::string name, std::vector<int> &ids, std::vector<double> &coefs, double &customExx);
};

} // namespace mrdft