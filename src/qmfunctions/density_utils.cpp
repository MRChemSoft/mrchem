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

#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/Parallel"
#include "MRCPP/MWFunctions"

#include "Density.h"
#include "Orbital.h"
#include "density_utils.h"
#include "orbital_utils.h"
#include <fstream>

using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/****************************************
 * Density related standalone functions *
 ****************************************/

namespace density {
Density compute(double prec, Orbital phi, DensityType spin);
void compute_local_X(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, DensityType spin);
void compute_local_XY(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, DensityType spin);
double compute_occupation(Orbital &phi, DensityType dens_spin);
} // namespace density

/** @brief Compute density as the square of an orbital
 *
 * This routine is similar to mrcpp::cplxfunc::multiply_real(), but it uses
 * mrcpp::square(phi) instead of mrcpp::multiply(phi, phi), which makes it
 * slightly faster.
 *
 */
Density density::compute(double prec, Orbital phi, DensityType spin) {
    double occ = density::compute_occupation(phi, spin);
    if (std::abs(occ) < mrcpp::MachineZero) return Density(false);

    Density rho(false);
    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        auto *real_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*real_2, phi.real());
        mrcpp::square(prec, *real_2, phi.real());
        sum_vec.push_back(std::make_tuple(occ, real_2));
    }
    if (phi.hasImag()) {
        auto *imag_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*imag_2, phi.imag());
        mrcpp::square(prec, *imag_2, phi.imag());
        sum_vec.push_back(std::make_tuple(occ, imag_2));
    }

    rho.alloc(NUMBER::Real);
    if (sum_vec.size() > 0) {
        mrcpp::build_grid(rho.real(), sum_vec);
        mrcpp::add(-1.0, rho.real(), sum_vec, 0);
        mrcpp::clear(sum_vec, true);
    } else {
        rho.real().setZero();
    }

    return rho;
}

/** @brief Compute density as the sum of squared orbitals
 *
 * MPI: Each rank first computes its own local density, which is then reduced
 *      to rank = 0 and broadcasted to all ranks.
 *
 */
void density::compute(double prec, Density &rho, OrbitalVector &Phi, DensityType spin) {
    int N_el = orbital::get_electron_number(Phi);
    double rel_prec = prec;        // prec for rho_i = |phi_i|^2
    double abs_prec = prec / N_el; // prec for rho = sum_i rho_i

    Density rho_loc(false);
    density::compute_local(rel_prec, rho_loc, Phi, spin);
    density::allreduce_density(abs_prec, rho, rho_loc);
}

/** @brief Compute transition density as rho = sum_i |x_i><phi_i| + |phi_i><y_i|
 *
 * MPI: Each rank first computes its own local density, which is then reduced
 *      to rank = 0 and broadcasted to all ranks. The rank distribution of Phi
 *      and X/Y must be the same.
 *
 */
void density::compute(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, DensityType spin) {
    int N_el = orbital::get_electron_number(Phi);
    double rel_prec = prec;        // prec for rho_i = |x_i><phi_i| + |phi_i><x_i|
    double abs_prec = prec / N_el; // prec for rho = sum_i rho_i

    // LUCA: rho_loc should get the "general" grid in order to make sure all densities are available on the same grid.
    Density rho_loc(false);
    if (&X == &Y) {
        density::compute_local_X(rel_prec, rho_loc, Phi, X, spin);
    } else {
        density::compute_local_XY(rel_prec, rho_loc, Phi, X, Y, spin);
    }
    density::allreduce_density(abs_prec, rho, rho_loc);
}

/** @brief Compute local density as the sum of own (MPI) orbitals
 */
void density::compute_local(double prec, Density &rho, OrbitalVector &Phi, DensityType spin) {
    int N_el = orbital::get_electron_number(Phi);
    double abs_prec = (mrcpp::mpi::numerically_exact) ? -1.0 : prec / N_el;
    if (not rho.hasReal()) rho.alloc(NUMBER::Real);

    if (rho.hasReal()) rho.real().setZero();
    if (rho.hasImag()) rho.imag().setZero();

    for (auto &phi_i : Phi) {
        if (mrcpp::mpi::my_orb(phi_i)) {
            Density rho_i = density::compute(prec, phi_i, spin);
            rho.add(1.0, rho_i); // Extends to union grid
            rho.crop(abs_prec);  // Truncates to given precision
        }
    }
}

/** @brief Compute local density as the sum of own (MPI) orbitals
 */
void density::compute_local(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, DensityType spin) {
    if (&X == &Y) {
        density::compute_local_X(prec, rho, Phi, X, spin);
    } else {
        density::compute_local_XY(prec, rho, Phi, X, Y, spin);
    }
}

void density::compute_local_X(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, DensityType spin) {
    int N_el = orbital::get_electron_number(Phi);
    double mult_prec = prec;       // prec for rho_i = |x_i><phi_i| + |phi_i><x_i|
    double add_prec = prec / N_el; // prec for rho = sum_i rho_i
    if (Phi.size() != X.size()) MSG_ERROR("Size mismatch");

    if (not rho.hasReal()) rho.alloc(NUMBER::Real);

    // Compute local density from own orbitals
    rho.real().setZero();
    for (int i = 0; i < Phi.size(); i++) {
        if (mrcpp::mpi::my_orb(Phi[i])) {
            if (not mrcpp::mpi::my_orb(X[i])) MSG_ABORT("Inconsistent MPI distribution");

            double occ = density::compute_occupation(Phi[i], spin);
            if (std::abs(occ) < mrcpp::MachineZero) continue; // next orbital if this one is not occupied!

            Density rho_i(false);
            mrcpp::cplxfunc::multiply_real(rho_i, Phi[i], X[i], mult_prec);
            rho.add(2.0 * occ, rho_i);
            rho.crop(add_prec);
        }
    }
}

void density::compute_local_XY(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, DensityType spin) {
    int N_el = orbital::get_electron_number(Phi);
    double mult_prec = prec;       // prec for rho_i = |x_i><phi_i| + |phi_i><y_i|
    double add_prec = prec / N_el; // prec for rho = sum_i rho_i
    if (Phi.size() != X.size()) MSG_ERROR("Size mismatch");
    if (Phi.size() != Y.size()) MSG_ERROR("Size mismatch");

    if (not rho.hasReal()) rho.alloc(NUMBER::Real);

    // Compute local density from own orbitals
    rho.real().setZero();
    for (int i = 0; i < Phi.size(); i++) {
        if (mrcpp::mpi::my_orb(Phi[i])) {
            if (not mrcpp::mpi::my_orb(X[i])) MSG_ABORT("Inconsistent MPI distribution");
            if (not mrcpp::mpi::my_orb(Y[i])) MSG_ABORT("Inconsistent MPI distribution");

            double occ = density::compute_occupation(Phi[i], spin);
            if (std::abs(occ) < mrcpp::MachineZero) continue; // next orbital if this one is not occupied!

            Density rho_x(false);
            Density rho_y(false);
            mrcpp::cplxfunc::multiply(rho_x, X[i], Phi[i].dagger(), mult_prec);
            mrcpp::cplxfunc::multiply(rho_y, Phi[i], Y[i].dagger(), mult_prec);
            rho.add(occ, rho_x);
            rho.add(occ, rho_y);
            rho.crop(add_prec);
        }
    }
}

void density::compute(double prec, Density &rho, mrcpp::GaussExp<3> &dens_exp) {
    if (not rho.hasReal()) rho.alloc(NUMBER::Real);
    mrcpp::build_grid(rho.real(), dens_exp);
    mrcpp::project(prec, rho.real(), dens_exp);
}

double density::compute_occupation(Orbital &phi, DensityType dens_spin) {
    double occ_a(0.0), occ_b(0.0), occ_p(0.0);
    if (phi.spin() == SPIN::Alpha) occ_a = (double)phi.occ();
    if (phi.spin() == SPIN::Beta) occ_b = (double)phi.occ();
    if (phi.spin() == SPIN::Paired) occ_p = (double)phi.occ();

    double occ(0.0);
    if (dens_spin == DensityType::Total) occ = occ_a + occ_b + occ_p;
    if (dens_spin == DensityType::Alpha) occ = occ_a + 0.5 * occ_p;
    if (dens_spin == DensityType::Beta) occ = occ_b + 0.5 * occ_p;
    if (dens_spin == DensityType::Spin) occ = occ_a - occ_b;

    return occ;
}

/** @brief Add up local density contributions and broadcast
 *
 * MPI: Each rank first computes its own local density, which is then reduced
 *      to rank = 0 and broadcasted to all ranks.
 *
 */
void density::allreduce_density(double prec, Density &rho_tot, Density &rho_loc) {
    // For numerically identical results in MPI we must first add
    // up orbital contributions onto their union grid, and THEN
    // crop the resulting density tree to the desired precision
    double part_prec = (mrcpp::mpi::numerically_exact) ? -1.0 : prec;
    // Add up local contributions into the grand master
    mrcpp::mpi::reduce_function(part_prec, rho_loc, mrcpp::mpi::comm_wrk);
    if (mrcpp::mpi::grand_master()) {
        // If numerically exact the grid is huge at this point
        if (mrcpp::mpi::numerically_exact) rho_loc.crop(prec);
    }

    if (not rho_tot.hasReal()) rho_tot.alloc(NUMBER::Real);

    if (rho_tot.isShared()) {
        int tag = 2002;
        if (mrcpp::mpi::share_master()) {
            // MPI grand master distributes to shared masters
            mrcpp::mpi::broadcast_function(rho_loc, mrcpp::mpi::comm_sh_group);
            // MPI shared masters copies the function into final memory
            mrcpp::copy_grid(rho_tot.real(), rho_loc.real());
            mrcpp::copy_func(rho_tot.real(), rho_loc.real());
        }
        // MPI share masters distributes to their sharing ranks
        mrcpp::mpi::share_function(rho_tot, 0, tag, mrcpp::mpi::comm_share);
    } else {
        // MPI grand master distributes to all ranks
        mrcpp::mpi::broadcast_function(rho_loc, mrcpp::mpi::comm_wrk);
       // All MPI ranks copies the function into final memory
        mrcpp::copy_grid(rho_tot.real(), rho_loc.real());
        mrcpp::copy_func(rho_tot.real(), rho_loc.real());
    }
}


// Function to read atomic density data from a file
void density::readAtomicDensity(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &rhoGrid) {
    std::vector<double> r, rho;
    std::ifstream file(path);
    std::string line;
    double r_, rho_;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> r_ >> rho_;
        if (rho_ > 0) {
            r.push_back(r_);
            rho.push_back(rho_);
        }
    }
    file.close();
    rGrid = Eigen::Map<Eigen::VectorXd>(r.data(), r.size());
    rhoGrid = Eigen::Map<Eigen::VectorXd>(rho.data(), rho.size());
}

} // namespace mrchem
