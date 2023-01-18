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

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "NuclearOperator.h"

#include "analyticfunctions/NuclearFunction.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/QMFunction.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/*! @brief NuclearOperator represents the function: sum_i Z_i/|r - R_i|
 *  @param nucs: Collection of nuclei that defines the potential
 *  @param proj_prec: Precision for projection of analytic function
 *  @param smooth_prec: Precision for smoothing of analytic function
 *  @param mpi_share: Should MPI ranks on the same machine share this function?
 */
NuclearOperator::NuclearOperator(const Nuclei &nucs, double proj_prec, double smooth_prec, bool mpi_share) {
    if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
    if (smooth_prec < 0.0) smooth_prec = proj_prec;

    Timer t_tot;
    mrcpp::print::header(2, "Projecting nuclear potential");

    // Setup local analytic function
    Timer t_loc;
    NuclearFunction f_loc;
    setupLocalPotential(f_loc, nucs, smooth_prec);

    // Scale precision by charge, since norm of potential is ~ to charge
    double Z_tot = 1.0 * chemistry::get_total_charge(nucs);
    double Z_loc = 1.0 * chemistry::get_total_charge(f_loc.getNuclei());
    double tot_prec = proj_prec / std::min(1.0 * Z_tot, std::sqrt(2.0 * Z_tot));
    double loc_prec = proj_prec / std::max(1.0, Z_loc); // relative prec

    // Scale precision by box size, so that accuracy is independent of box size
    double vol = 1.0;
    for (int i = 0; i < 3; i++) vol *= MRA->getWorldBox().getBoxLength(i);
    vol /= 262144;                   // we use as reference a cube 64x64x64
    vol = std::max(1.0, vol);        // do not scale for smaller boxes
    loc_prec /= pow(vol, 1.0 / 6.0); // norm of 1/r over the box ~ root_6(Volume)

    // Project local potential
    QMFunction V_loc(false);
    qmfunction::project(V_loc, f_loc, NUMBER::Real, loc_prec);
    t_loc.stop();
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local potential", V_loc, t_loc);

    // Collect local potentials
    Timer t_com;
    auto V_tot = std::make_shared<QMPotential>(1, mpi_share);
    allreducePotential(tot_prec, *V_tot, V_loc);
    t_com.stop();

    t_tot.stop();
    print_utils::qmfunction(2, "Allreduce potential", *V_tot, t_com);
    mrcpp::print::footer(2, t_tot, 2);

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_tot;
    O.name() = "V_nuc";
}

void NuclearOperator::setupLocalPotential(NuclearFunction &f_loc, const Nuclei &nucs, double smooth_prec) const {
    int pprec = Printer::getPrecision();
    int w0 = Printer::getWidth() - 1;
    int w1 = 5;
    int w2 = 8;
    int w3 = 2 * w0 / 9;
    int w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::string(w4, ' ');
    o_head << std::setw(w3) << "Charge";
    o_head << std::setw(w3) << "Precision";
    o_head << std::setw(w3) << "Smoothing";

    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double c = detail::nuclear_potential_smoothing(smooth_prec, Z);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mpi::numerically_exact) ? 0 : k % mpi::orb_size;
        if (mpi::orb_rank == proj_rank) f_loc.push_back(nuc, c);

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << nuc.getElement().getSymbol();
        o_row << std::string(w4, ' ');
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << Z;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth_prec;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << c;
        println(2, o_row.str());
    }
}

void NuclearOperator::allreducePotential(double prec, QMFunction &V_tot, QMFunction &V_loc) const {
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
