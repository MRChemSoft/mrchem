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

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "SmearedNuclearOperator.h"

#include "analyticfunctions/NuclearFunction.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/QMFunction.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "utils/math_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/*! @brief SmearedNuclearOperator represents the function: sum_i Z_i/|r - R_i|
 *  @param nucs: Collection of nuclei that defines the potential
 *  @param proj_prec: Precision for projection of analytic function
 *  @param smooth_prec: Precision for smoothing of analytic function
 *  @param mpi_share: Should MPI ranks on the same machine share this function?
 */
SmearedNuclearOperator::SmearedNuclearOperator(Nuclei nucs, double proj_prec, double smooth_prec, double rc, bool mpi_share, std::shared_ptr<OrbitalVector> Phi) {
    if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
    if (smooth_prec < 0.0) smooth_prec = proj_prec;

    QMFunction V_loc(false);
    setupLocalPotential(V_loc, nucs, proj_prec, smooth_prec, rc);

    // Scale precision by system size
    double Z_tot = 1.0 * chemistry::get_total_charge(nucs);
    double tot_prec = proj_prec / std::min(1.0 * Z_tot, std::sqrt(2.0 * Z_tot));
    // Project local potential
    auto V_tot = std::make_shared<QMPotential>(1, mpi_share);
    allreducePotential(tot_prec, *V_tot, V_loc);

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_tot;
}

void SmearedNuclearOperator::setupLocalPotential(QMFunction &V_loc, const Nuclei &nucs, double proj_prec, double smooth_prec, double rc) const {

    NuclearFunction loc_func;
    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double c = detail::nuclear_potential_smoothing(smooth_prec, Z);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mpi::numerically_exact) ? 0 : k % mpi::orb_size;
        if (mpi::orb_rank == proj_rank) loc_func.push_back(nuc, c);
    }

    auto V_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
        auto v_g = 0.0;
        for (auto &nuc : nucs) {
            auto R = math_utils::calc_distance(r, nuc.getCoord());
            auto v_i = 0.0;
            if (R <= rc and R >= 0) {
                v_i = -(9.0 * std::pow(R, 7.0) - 30.0 * std::pow(R, 6.0) * rc + 28.0 * std::pow(R, 5.0) * rc * rc - 14.0 * R * R * std::pow(rc, 5.0) + 12.0 * std::pow(rc, 7.0)) /
                      (5.0 * std::pow(rc, 8.0));
            } else {
                v_i = -1.0 / R;
            }
            v_g += v_i * nuc.getCharge();
        }

        return v_g;
    };
    // Scale precision by system size
    int Z_tot = chemistry::get_total_charge(nucs);
    double abs_prec = proj_prec / Z_tot;

    QMFunction V_nuc_loc(false);
    QMFunction V_smear_loc(false);

    qmfunction::project(V_nuc_loc, loc_func, NUMBER::Real, abs_prec);
    qmfunction::project(V_smear_loc, V_smear, NUMBER::Real, abs_prec);

    if (not V_loc.hasReal()) V_loc.alloc(NUMBER::Real);
    qmfunction::add(V_loc, 1.0, V_nuc_loc, -1.0, V_smear_loc, -1);
}

void SmearedNuclearOperator::allreducePotential(double prec, QMFunction &V_tot, QMFunction &V_loc) const {
    // Add up local contributions into the grand master
    mpi::reduce_function(prec, V_loc, mpi::comm_orb);
    if (mpi::grand_master()) {
        // If numerically exact the grid is huge at this point
        if (mpi::numerically_exact) V_loc.crop(prec);
    }

    if (not V_tot.hasReal()) V_tot.alloc(NUMBER::Real);
    if (V_tot.isShared()) {
        int tag = 3141;
        // MPI grand master distributes to shared masters
        mpi::broadcast_function(V_loc, mpi::comm_sh_group);
        if (mpi::share_master()) {
            // MPI shared masters copies the function into final memory
            mrcpp::copy_grid(V_tot.real(), V_loc.real());
            mrcpp::copy_func(V_tot.real(), V_loc.real());
        }
        // MPI share masters distributes to their sharing ranks
        mpi::share_function(V_tot, 0, tag, mpi::comm_share);
    } else {
        // MPI grand master distributes to all ranks
        mpi::broadcast_function(V_loc, mpi::comm_orb);
        // All MPI ranks copies the function into final memory
        mrcpp::copy_grid(V_tot.real(), V_loc.real());
        mrcpp::copy_func(V_tot.real(), V_loc.real());
    }
}

} // namespace mrchem
