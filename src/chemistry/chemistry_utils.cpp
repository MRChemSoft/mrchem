/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "chemistry_utils.h"
#include "Nucleus.h"
#include "PhysicalConstants.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/density_utils.h"
#include "qmoperators/one_electron/H_E_dip.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"

#include <MRCPP/Gaussians>

namespace mrchem {

/** @brief computes the repulsion self energy of a set of nuclei
 *
 * @param[in] nucs the set of nuclei
 *
 */
double chemistry::compute_nuclear_repulsion(const Nuclei &nucs) {
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int i = 0; i < nNucs; i++) {
        const Nucleus &nuc_i = nucs[i];
        const double Z_i = nuc_i.getCharge();
        const mrcpp::Coord<3> &R_i = nuc_i.getCoord();
        for (int j = i + 1; j < nNucs; j++) {
            const Nucleus &nuc_j = nucs[j];
            const double Z_j = nuc_j.getCharge();
            const mrcpp::Coord<3> &R_j = nuc_j.getCoord();
            double R_ij = math_utils::calc_distance(R_i, R_j);
            E_nuc += (Z_i * Z_j) / R_ij;
        }
    }
    return E_nuc;
}

/** @brief Returns the sum of atomic charges*/
double chemistry::get_total_charge(const Nuclei &nucs) {
    double charge = 0;
    for (const auto &nuc : nucs) charge += nuc.getCharge();
    return charge;
}

/** @breif computes the nuclear density as a sum of narrow Gaussians */
Density chemistry::compute_nuclear_density(double prec, const Nuclei &nucs, double alpha) {
    auto beta = std::pow(alpha / mrcpp::pi, 3.0 / 2.0);
    auto gauss = mrcpp::GaussExp<3>();
    for (auto i = 0; i < nucs.size(); i++) {
        const auto &nuc_i = nucs[i];
        const auto Z_i = nuc_i.getCharge();
        const auto &R_i = nuc_i.getCoord();
        auto gauss_f = mrcpp::GaussFunc<3>(alpha, beta * Z_i, R_i);
        gauss.append(gauss_f);
    }

    if ((*MRA).getWorldBox().isPeriodic()) {
        auto period = (*MRA).getWorldBox().getScalingFactors();
        gauss.periodify(period);
    }
    Density rho(false);
    density::compute(prec, rho, gauss);
    return rho;
}
Density chemistry::compute_nuclear_density_smeared(double prec, Nuclei nucs, double rc, double period) {
    Density rho(false);
    if (not rho.hasReal()) rho.alloc(NUMBER::Real);

    nucs = periodic::periodify_nuclei(nucs, period);

    auto b_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
        auto g_rc = 0.0;
        for (auto &nuc : nucs) {
            auto R = math_utils::calc_distance(r, nuc.getCoord());
            auto g_i = 0.0;
            if (R <= rc and R >= 0) { g_i = -21.0 * std::pow((R - rc), 3.0) * (6.0 * R * R + 3.0 * R * rc + rc * rc) / (5.0 * mrcpp::pi * std::pow(rc, 8.0)); }
            g_rc += g_i * nuc.getCharge();
        }
        return g_rc;
    };

    // We make fictisious gaussians to build a grid, this makes our project
    // find b_smeared.
    auto beta = 1.0e4;
    auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
    mrcpp::GaussExp<3> f_exp;
    for (auto &nuc : nucs) {
        auto g_func = mrcpp::GaussFunc<3>(beta, alpha, nuc.getCoord());
        f_exp.append(g_func);
    }

    mrcpp::build_grid(rho.real(), f_exp);
    mrcpp::project<3>(prec, rho.real(), b_smear);
    return rho;
}
/** @breif computes the nuclear density as a sum of narrow Gaussians */
double chemistry::compute_nuclear_self_repulsion(const Nuclei &nucs, double alpha) {
    auto beta = std::pow(alpha / MATHCONST::pi, 3.0 / 2.0);

    auto self_rep = 0.0;
    for (auto i = 0; i < nucs.size(); i++) {
        const auto &nuc_i = nucs[i];
        const auto Z_i = nuc_i.getCharge();
        const auto &R_i = nuc_i.getCoord();
        auto gauss_f = mrcpp::GaussFunc<3>(alpha, beta, R_i);
        self_rep += Z_i * Z_i * gauss_f.calcCoulombEnergy(gauss_f);
    }
    return self_rep;
}

} // namespace mrchem
