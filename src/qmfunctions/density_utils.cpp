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

void density::compute(double prec, Density &rho, Orbital phi, int spin) {
    double occ_a(0.0), occ_b(0.0), occ_p(0.0);
    if (phi.spin() == SPIN::Alpha)  occ_a = (double) phi.occ();
    if (phi.spin() == SPIN::Beta)   occ_b = (double) phi.occ();
    if (phi.spin() == SPIN::Paired) occ_p = (double) phi.occ();

    double occ(0.0);
    if (spin == DENSITY::Total) occ = occ_a + occ_b + occ_p;
    if (spin == DENSITY::Alpha) occ = occ_a + 0.5*occ_p;
    if (spin == DENSITY::Beta)  occ = occ_b + 0.5*occ_p;
    if (spin == DENSITY::Spin)  occ = occ_a - occ_b;

    if (std::abs(occ) < mrcpp::MachineZero) {
        rho.real().setZero();
        return;
    }

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*real_2, phi.real());
        mrcpp::square(prec, *real_2, phi.real());
        sum_vec.push_back(std::make_tuple(occ, real_2));
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*imag_2, phi.imag());
        mrcpp::square(prec, *imag_2, phi.imag());
        sum_vec.push_back(std::make_tuple(occ, imag_2));
    }
    mrcpp::build_grid(rho.real(), sum_vec);
    mrcpp::add(-1.0, rho.real(), sum_vec, 0);
    mrcpp::clear(sum_vec, true);
}

void density::compute(double prec, Density &rho, OrbitalVector &Phi, int spin) {
    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i

    FunctionTreeVector<3> dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            Density rho_i;
            rho_i.alloc(NUMBER::Real);
            mrcpp::copy_grid(rho_i.real(), rho.real());
            density::compute(mult_prec, rho_i, Phi[i], spin);
            dens_vec.push_back(std::make_tuple(1.0, &(rho_i.real())));
            rho_i.clear(); // release FunctionTree pointers to dens_vec
        }
    }

    if (add_prec > 0.0) {
        mrcpp::add(add_prec, rho.real(), dens_vec);
    } else if (dens_vec.size() > 0) {
        mrcpp::build_grid(rho.real(), dens_vec);
        mrcpp::add(-1.0, rho.real(), dens_vec, 0);
    }
    mrcpp::clear(dens_vec, true);

    mpi::reduce_density(rho, mpi::comm_orb);
    mpi::broadcast_density(rho, mpi::comm_orb);
}

void density::compute(double prec, Density &rho, mrcpp::GaussExp<3> &dens_exp, int spin) {
    rho.alloc(NUMBER::Real);
    rho.setSpin(spin);
    mrcpp::project(prec, rho.real(), dens_exp);
}

} //namespace mrchem
