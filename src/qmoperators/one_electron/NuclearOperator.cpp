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

#include <MRCPP/MWOperators>
#include <MRCPP/Parallel>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "NuclearOperator.h"

#include<fstream>
#include<limits>
#include<nlohmann/json.hpp>

#include "analyticfunctions/NuclearFunction.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/QMFunction.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "utils/print_utils.h"
#include "utils/math_utils.h"


using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {


NuclearOperator::NuclearOperator(const Nuclei &nucs, double proj_prec, bool mpi_share, int nuc_model) {
    if (nuc_model == 1) {
        projectHFYGB(nucs, proj_prec, mpi_share);
    } else if (nuc_model == 2) {
        projectHomogeneusSphere(nucs, proj_prec, mpi_share);
    } else if (nuc_model == 3) {
        projectGaussian(nucs, proj_prec, mpi_share);
    } else {
        NOT_REACHED_ABORT("Invalid nuclear type");
    }
}

/** Compute finite nucleus potential by projecting analytic expression for the potential.
*/
void NuclearOperator::projectHFYGB(const Nuclei &nucs, double proj_prec, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(0, "Projecting HFYGB Nuclear Potential");
 
    auto f_smooth = [&nucs,proj_prec] (const mrcpp::Coord<3> &r) -> double {
        double c = -1.0 / (3.0 * mrcpp::root_pi);
        double tmp_sum = 0.0;
        for (auto &nuc : nucs) {
            const auto Z_I = nuc.getCharge();
            const auto &R_I = nuc.getCoord();
            const double Z_I_2 = Z_I * Z_I;
            const double Z_I_5 = Z_I_2 * Z_I_2 * Z_I;
            const double factor = std::cbrt(0.00435 * proj_prec / Z_I_5);
            const double R = math_utils::calc_distance(r, R_I) / factor;
            const double u = -std::erf(R)/R + c*(std::exp(-(R*R)) + 16.0*std::exp(-4.0*R*R));
            tmp_sum += (Z_I * u) / factor;
        }
        return tmp_sum;
    };

    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    qmfunction::project(*V_nuc, f_smooth, NUMBER::Real, proj_prec);
    t_pot.stop();

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}

void NuclearOperator::projectHomogeneusSphere(const Nuclei &nucs, double proj_prec, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(0, "Projecting Homogeneous Sphere Nuclear Potential");

    std::ifstream f("param_V.json");
    while (f.peek() == '#') f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    nlohmann::json param_V = nlohmann::json::parse(f);
    f.close();
 
    auto f_sphere = [&nucs,&param_V] (const mrcpp::Coord<3> &r) -> double {
        double tmp_sum = 0.0;
        for (auto &nuc : nucs) {
            const auto &A = nuc.getElement().getSymbol();
            double RMS = 0.0;
            for (auto item : param_V["element"]) {
                const auto name = item["name"];
                if (name == A) {
                    RMS = item["rms"];
                    break;
                }
            }
            const auto RMS2 = RMS*RMS;
            const auto Z_I = nuc.getCharge();
            const auto &R_I = nuc.getCoord();
            double R = math_utils::calc_distance(r, R_I);
            double R0 = std::sqrt(RMS2*(5.0/3.0));
            double prec, factor;
            if (R <= R0) {
                prec = -Z_I / (2.0*R0);
                factor = 3.0 - (R*R)/(R0*R0);
            } else { 
                prec = -Z_I / R;
                factor = 1.0;
            }
            tmp_sum += prec * factor;
        }
        return tmp_sum;
    };

    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    qmfunction::project(*V_nuc, f_sphere, NUMBER::Real, proj_prec);
    t_pot.stop();

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}

void NuclearOperator::projectGaussian(const Nuclei &nucs, double proj_prec, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(0, "Projecting Gaussian Nuclear Potential");

    std::ifstream f("param_V.json");
    while (f.peek() == '#') f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    nlohmann::json param_V = nlohmann::json::parse(f);
    f.close();
 
    auto f_smooth = [&nucs,&param_V] (const mrcpp::Coord<3> &r) -> double {
        double tmp_sum = 0.0;
        for (auto &nuc : nucs) {
            const auto &A = nuc.getElement().getSymbol();
            double epsilon = 0.0;
            for (auto item : param_V["element"]) {
                const auto name = item["name"];
                if (name == A) {
                    epsilon = item["xi"];
                    break;
                }
            }
            const auto Z_I = nuc.getCharge();
            const auto &R_I = nuc.getCoord();
            const double R = math_utils::calc_distance(r, R_I);
            const double prec = -Z_I / R;
            const double u = std::erf(std::sqrt(epsilon) * R);
            tmp_sum += prec * u;
        }
        return tmp_sum;
    };

    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    qmfunction::project(*V_nuc, f_smooth, NUMBER::Real, proj_prec);
    t_pot.stop();

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}

/** Compute finite nucleus potential by projecting Gaussian charge distribution and then apply Poisson.
*   Should work for any number of nuclei, but all get the same exponent.
*/
void NuclearOperator::applyGaussian(const Nuclei &nucs, double proj_prec, double apply_prec, double exponent, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(2, "Projecting nuclear density");
    println(0, "Nuclear exponent: " << exponent);

    mrcpp::PoissonOperator P(*MRA, proj_prec);

    Timer t_rho;
    auto rho_nuc = chemistry::compute_nuclear_density(proj_prec, nucs, exponent);
    rho_nuc.rescale(-1.0);
    t_rho.stop();

    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    V_nuc->alloc(NUMBER::Real);
    mrcpp::apply(apply_prec, V_nuc->real(), P, rho_nuc.real());
    t_pot.stop();

    t_tot.stop();
    print_utils::qmfunction(2, "Apply Poisson", *V_nuc, t_pot);
    mrcpp::print::footer(2, t_tot, 2);

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}

/*! @brief NuclearOperator represents the function: sum_i Z_i/|r - R_i|
 *  @param nucs: Collection of nuclei that defines the potential
 *  @param proj_prec: Precision for projection of analytic function
 *  @param smooth_prec: Precision for smoothing of analytic function
 *  @param mpi_share: Should MPI ranks on the same machine share this function?
 */
NuclearOperator::NuclearOperator(const Nuclei &nucs, double proj_prec, double smooth_prec, bool mpi_share) {
    if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
    if (smooth_prec < 0.0) smooth_prec = proj_prec;
    this->nucs = nucs;

    Timer t_tot;
    mrcpp::print::header(2, "Projecting nuclear potential");

    // Setup local analytic function
    Timer t_loc;
    NuclearFunction f_loc;
    setupLocalPotential(f_loc, nucs, smooth_prec, false);
    Nuc_func = NuclearFunction(nucs, smooth_prec);

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
    mrcpp::ComplexFunction V_loc(false);
    mrcpp::cplxfunc::project(V_loc, f_loc, NUMBER::Real, loc_prec);
    t_loc.stop();
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local potential", V_loc, t_loc);

    // Collect local potentials
    Timer t_com;

    auto V_tot = std::make_shared<QMPotential>(1, mpi_share);
    allreducePotential(tot_prec, *V_tot, V_loc);
    V_func = *V_tot;

    t_com.stop();

    t_tot.stop();
    print_utils::qmfunction(2, "Allreduce potential", *V_tot, t_com);
    mrcpp::print::footer(2, t_tot, 2);

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_tot;
    O.name() = "V_nuc";
}

void NuclearOperator::setupLocalPotential(NuclearFunction &f_loc, const Nuclei &nucs, double smooth_prec, bool print) const {
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

    if (print) println(2, o_head.str());
    if (print) mrcpp::print::separator(2, '-');

    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double c = detail::nuclear_potential_smoothing(smooth_prec, Z);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mrcpp::mpi::numerically_exact) ? 0 : k % mrcpp::mpi::wrk_size;
        if (mrcpp::mpi::wrk_rank == proj_rank) f_loc.push_back(nuc, c);

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << nuc.getElement().getSymbol();
        o_row << std::string(w4, ' ');
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << Z;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth_prec;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << c;
        if (print) println(2, o_row.str());
    }
}

void NuclearOperator::allreducePotential(double prec, mrcpp::ComplexFunction &V_tot, mrcpp::ComplexFunction &V_loc) const {
    // Add up local contributions into the grand master
    mrcpp::mpi::reduce_function(prec, V_loc, mrcpp::mpi::comm_wrk);
    if (mrcpp::mpi::grand_master()) {
        // If numerically exact the grid is huge at this point
        if (mrcpp::mpi::numerically_exact) V_loc.crop(prec);
    }

    if (not V_tot.hasReal()) V_tot.alloc(NUMBER::Real);
    if (V_tot.isShared()) {
        int tag = 3141;
        // MPI grand master distributes to shared masters
        mrcpp::mpi::broadcast_function(V_loc, mrcpp::mpi::comm_sh_group);
        if (mrcpp::mpi::share_master()) {
            // MPI shared masters copies the function into final memory
            mrcpp::copy_grid(V_tot.real(), V_loc.real());
            mrcpp::copy_func(V_tot.real(), V_loc.real());
        }
        // MPI share masters distributes to their sharing ranks
        mrcpp::mpi::share_function(V_tot, 0, tag, mrcpp::mpi::comm_share);
    } else {
        // MPI grand master distributes to all ranks
        mrcpp::mpi::broadcast_function(V_loc, mrcpp::mpi::comm_wrk);
        // All MPI ranks copies the function into final memory
        mrcpp::copy_grid(V_tot.real(), V_loc.real());
        mrcpp::copy_func(V_tot.real(), V_loc.real());
    }
}

} // namespace mrchem
