/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/Gaussians"

#include "parallel.h"

#include "QMFunction.h"
#include "qmfunction_utils.h"
#include "Density.h"
#include "density_utils.h"
#include "Orbital.h"
#include "orbital_utils.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/****************************************
 * Density related standalone functions *
 ****************************************/

namespace density {
void compute(double prec, Density &rho, Orbital phi, int spin);
double compute_occupation(Orbital &phi, int dens_spin);
}

/** @brief Compute density as the square of an orbital
 *
 * This routine is similar to qmfunction::multiply_real(), but it uses
 * mrcpp::square(phi) instead of mrcpp::multiply(phi, phi), which makes it
 * slightly faster.
 *
 */
void density::compute(double prec, Density &rho, Orbital phi, int spin) {
    double occ = density::compute_occupation(phi, spin);
    bool phi_contributes = (std::abs(occ) > mrcpp::MachineZero);

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal() and phi_contributes) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*real_2, phi.real());
        mrcpp::square(prec, *real_2, phi.real());
        sum_vec.push_back(std::make_tuple(occ, real_2));
    }
    if (phi.hasImag() and phi_contributes) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*imag_2, phi.imag());
        mrcpp::square(prec, *imag_2, phi.imag());
        sum_vec.push_back(std::make_tuple(occ, imag_2));
    }

    if (sum_vec.size() > 0) {
        if (not rho.hasReal()) rho.alloc(NUMBER::Real);
        mrcpp::build_grid(rho.real(), sum_vec);
        mrcpp::add(-1.0, rho.real(), sum_vec, 0);
        mrcpp::clear(sum_vec, true);
    } else if (rho.hasReal()) {
        rho.real().setZero();
    }
    if (rho.hasImag()) {
        rho.imag().setZero();
    }
}

/** @brief Compute density as the sum of squared orbitals
 *
 * MPI: Each rank first computes its own local density, which is then reduced
 *      to rank = 0 and broadcasted to all ranks.
 *
 */
void density::compute(double prec, Density &rho, OrbitalVector &Phi, int spin) {
    int N_el = orbital::get_electron_number(Phi);
    double mult_prec = prec;     // prec for rho_i = |phi_i|^2
    double add_prec = prec/N_el; // prec for rho = sum_i rho_i

    FunctionTreeVector<3> dens_vec;
    for (auto &phi_i : Phi) {
        if (mpi::my_orb(phi_i)) {
            Density rho_i;
            density::compute(mult_prec, rho_i, phi_i, spin);
            if (rho_i.hasReal()) dens_vec.push_back(std::make_tuple(1.0, &(rho_i.real())));
            if (rho_i.hasImag()) MSG_ERROR("Density should be real");
        }
    }

    if (dens_vec.size() > 0) {
        if (not rho.hasReal()) rho.alloc(NUMBER::Real);
        mrcpp::add(add_prec, rho.real(), dens_vec);
    } else if (rho.hasReal()) {
        rho.real().setZero();
    }
    if (rho.hasImag()) {
        rho.imag().setZero();
    }

    mrcpp::clear(dens_vec, true);

    mpi::reduce_density(rho, mpi::comm_orb);
    mpi::broadcast_density(rho, mpi::comm_orb);
}

/** @brief Compute transition density as rho = sum_i |x_i><phi_i| + |phi_i><x_i|
 *
 * Exploits the fact that the resulting density must be real.
 *
 * MPI: Each rank first computes its own local density, which is then reduced
 *      to rank = 0 and broadcasted to all ranks. The rank distribution of Phi
 *      and X must be the same.
 *
 */
void density::compute(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, int spin) {
    int N_el = orbital::get_electron_number(Phi);
    double mult_prec = prec;     // prec for rho_i = |x_i><phi_i| + |phi_i><x_i|
    double add_prec = prec/N_el; // prec for rho = sum_i rho_i
    if (Phi.size() != X.size()) MSG_ERROR("Size mismatch");

    FunctionTreeVector<3> dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            if (not mpi::my_orb(X[i])) MSG_FATAL("Inconsistent MPI distribution");

            double occ = density::compute_occupation(Phi[i], spin);
            if (std::abs(occ) < mrcpp::MachineZero) continue; //next orbital if this one is not occupied!

            Density rho_i;
            qmfunction::multiply_real(rho_i, Phi[i], X[i], mult_prec);
            if (rho_i.hasReal()) dens_vec.push_back(std::make_tuple(2.0*occ, &(rho_i.real())));
        }
    }

    if (dens_vec.size() > 0) {
        if (not rho.hasReal()) rho.alloc(NUMBER::Real);
        mrcpp::add(add_prec, rho.real(), dens_vec);
    } else if (rho.hasReal()) {
        rho.real().setZero();
    }
    if (rho.hasImag()) {
        rho.imag().setZero();
    }

    mrcpp::clear(dens_vec, true);

    mpi::reduce_density(rho, mpi::comm_orb);
    mpi::broadcast_density(rho, mpi::comm_orb);
}

/** @brief Compute transition density as rho = sum_i |x_i><phi_i| + |phi_i><y_i|
 *
 * MPI: Each rank first computes its own local density, which is then reduced
 *      to rank = 0 and broadcasted to all ranks. The rank distribution of Phi
 *      and X/Y must be the same.
 *
 */
void density::compute(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, int spin) {
    int N_el = orbital::get_electron_number(Phi);
    double mult_prec = prec;     // prec for rho_i = |x_i><phi_i| + |phi_i><y_i|
    double add_prec = prec/N_el; // prec for rho = sum_i rho_i
    if (Phi.size() != X.size()) MSG_ERROR("Size mismatch");
    if (Phi.size() != Y.size()) MSG_ERROR("Size mismatch");

    FunctionTreeVector<3> dens_real;
    FunctionTreeVector<3> dens_imag;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            if (not mpi::my_orb(X[i])) MSG_FATAL("Inconsistent MPI distribution");
            if (not mpi::my_orb(Y[i])) MSG_FATAL("Inconsistent MPI distribution");

            double occ = density::compute_occupation(Phi[i], spin);
            if (std::abs(occ) < mrcpp::MachineZero) continue; //next orbital if this one is not occupied!

            Density rho_x;
            Density rho_y;
            qmfunction::multiply(rho_x, X[i], Phi[i].dagger(), mult_prec);
            qmfunction::multiply(rho_y, Phi[i], Y[i].dagger(), mult_prec);

            if (rho_x.hasReal()) dens_real.push_back(std::make_tuple(occ, &(rho_x.real())));
            if (rho_y.hasReal()) dens_real.push_back(std::make_tuple(occ, &(rho_y.real())));
            if (rho_x.hasImag()) dens_imag.push_back(std::make_tuple(occ, &(rho_x.imag())));
            if (rho_y.hasImag()) dens_imag.push_back(std::make_tuple(occ, &(rho_y.imag())));
        }
    }

    if (dens_real.size() > 0) {
        if (not rho.hasReal()) rho.alloc(NUMBER::Real);
        mrcpp::add(add_prec,  rho.real(), dens_real);
    } else if (rho.hasReal()) {
        rho.real().setZero();
    }

    if (dens_imag.size() > 0) {
        if (not rho.hasImag()) rho.alloc(NUMBER::Imag);
        mrcpp::add(add_prec,  rho.imag(), dens_imag);
    } else if (rho.hasImag()) {
        rho.imag().setZero();
    }

    mrcpp::clear(dens_real, true);
    mrcpp::clear(dens_imag, true);

    mpi::reduce_density(rho, mpi::comm_orb);
    mpi::broadcast_density(rho, mpi::comm_orb);
}

void density::compute(double prec, Density &rho, mrcpp::GaussExp<3> &dens_exp, int spin) {
    rho.alloc(NUMBER::Real);
    mrcpp::project(prec, rho.real(), dens_exp);
}

double density::compute_occupation(Orbital &phi, int dens_spin) {
    double occ_a(0.0), occ_b(0.0), occ_p(0.0);
    if (phi.spin() == SPIN::Alpha)  occ_a = (double) phi.occ();
    if (phi.spin() == SPIN::Beta)   occ_b = (double) phi.occ();
    if (phi.spin() == SPIN::Paired) occ_p = (double) phi.occ();
    
    double occ(0.0);
    if (dens_spin == DENSITY::Total) occ = occ_a + occ_b + occ_p;
    if (dens_spin == DENSITY::Alpha) occ = occ_a + 0.5*occ_p;
    if (dens_spin == DENSITY::Beta)  occ = occ_b + 0.5*occ_p;
    if (dens_spin == DENSITY::Spin)  occ = occ_a - occ_b;

    return occ;
}
    
} //namespace mrchem
