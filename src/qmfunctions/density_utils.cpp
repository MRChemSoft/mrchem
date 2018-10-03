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

void compute(double prec, Density &rho, Orbital phi, double occ);
void compute(double prec, Density &rho, Orbital ket, Orbital bra, double coeff, int type);
double compute_occupation(int orb_spin, int orb_occ, int dens_spin);

void density::compute(double prec, Density &rho, Orbital phi, double occ) {

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

void density::compute(double prec, Density &rho, Orbital ket, Orbital bra, double coeff, int type) {

    double ket_conj(1.0), bra_conj(1.0);
    if (ket.conjugate()) ket_conj = -1.0;
    if (bra.conjugate()) bra_conj = -1.0;
    
    if (type == NUMBER::Total) {
        qmfunction::multiply(ket, ket_conj, bra, -1.0 * bra_conj, rho, prec);
        rho.real().rescale(coeff);
        rho.imag().rescale(coeff);
    } else if (type == NUMBER::Real) {
        qmfunction::multiply_real(ket, ket_conj, bra, -1.0 * bra_conj, rho, prec);
        rho.real().rescale(coeff);
    } else {
        MSG_FATAL("No such case");
    }
}

void density::compute(double prec, Density &rho, OrbitalVector &Phi, int spin) {
    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i

    FunctionTreeVector<3> dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            double occ = compute_occupation(Phi[i].spin(), Phi[i].occ(), spin);
            if (std::abs(occ) < mrcpp::MachineZero) continue;
            Density rho_i;
            rho_i.alloc(NUMBER::Real);
            mrcpp::copy_grid(rho_i.real(), rho.real());
            density::compute(mult_prec, rho_i, Phi[i], occ);
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


//LUCA Is the MPI distribution of Phi identical to Phi_x and Phi_y??? Now I 
void density::compute(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &Phi_x, int spin) {
    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i
    if (Phi.size() != Phi_x.size()) MSG_ERROR("Size mismatch");
    
    FunctionTreeVector<3> dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            double occ = compute_occupation(Phi[i].spin(), Phi[i].occ(), spin);
            if (std::abs(occ) < mrcpp::MachineZero) continue; //next orbital if this one is not occupied!
            Density *rho_i = new Density(); 
            rho_i->alloc(NUMBER::Real);
            mrcpp::copy_grid(rho_i->real(), rho.real());
            density::compute(mult_prec, *rho_i, Phi[i], Phi_x[i], 2.0 * occ, NUMBER::Real);
            dens_vec.push_back(std::make_tuple(1.0, &(rho_i->real())));
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

//LUCA Is the MPI distribution of Phi identical to Phi_x and Phi_y??? Now I 
void density::compute(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &Phi_x, OrbitalVector &Phi_y, int spin) {
    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i
    if (Phi.size() != Phi_x.size()) MSG_ERROR("Size mismatch");
    
    FunctionTreeVector<3> dens_real;
    FunctionTreeVector<3> dens_imag;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            //if(mpi::my_orb(Phi_x[i])) ....
            double occ = compute_occupation(Phi[i].spin(), Phi[i].occ(), spin);
            if (std::abs(occ) < mrcpp::MachineZero) continue; //next orbital if this one is not occupied!
            Density *rho_x = new Density(); 
            Density *rho_y = new Density(); 
            rho_x->alloc(NUMBER::Real);
            rho_y->alloc(NUMBER::Real);
            rho_x->alloc(NUMBER::Imag);
            rho_y->alloc(NUMBER::Imag);
            mrcpp::copy_grid(rho_x->real(), rho.real());
            mrcpp::copy_grid(rho_y->real(), rho.real());
            mrcpp::copy_grid(rho_x->imag(), rho.imag());
            mrcpp::copy_grid(rho_y->imag(), rho.imag());
            density::compute(mult_prec, *rho_x, Phi_x[i], Phi[i], occ, NUMBER::Total);
            density::compute(mult_prec, *rho_y, Phi[i], Phi_y[i], occ, NUMBER::Total);
            dens_real.push_back(std::make_tuple(1.0, &(rho_x->real())));
            dens_real.push_back(std::make_tuple(1.0, &(rho_y->real())));
            dens_imag.push_back(std::make_tuple(1.0, &(rho_x->imag())));
            dens_imag.push_back(std::make_tuple(1.0, &(rho_y->imag())));
        }
    }

    if (add_prec > 0.0) {
        mrcpp::add(add_prec, rho.real(), dens_real);
        mrcpp::add(add_prec, rho.imag(), dens_imag);
    } else {
        if (dens_real.size() > 0) {
            mrcpp::build_grid(rho.real(), dens_real);
            mrcpp::add(-1.0,  rho.real(), dens_real, 0);
        }
        if (dens_imag.size() > 0) {
            mrcpp::build_grid(rho.imag(), dens_imag);
            mrcpp::add(-1.0,  rho.imag(), dens_imag, 0);
        }
    }
    mrcpp::clear(dens_real, true);
    mrcpp::clear(dens_imag, true);

    mpi::reduce_density(rho, mpi::comm_orb);
    mpi::broadcast_density(rho, mpi::comm_orb);
}

void density::compute(double prec, Density &rho, mrcpp::GaussExp<3> &dens_exp, int spin) {
    rho.alloc(NUMBER::Real);
    rho.setSpin(spin);
    mrcpp::project(prec, rho.real(), dens_exp);
}

double density::compute_occupation(int orb_spin, int orb_occ, int dens_spin) {
    double occ_a(0.0), occ_b(0.0), occ_p(0.0);
    if (orb_spin == SPIN::Alpha)  occ_a = (double) orb_occ;
    if (orb_spin == SPIN::Beta)   occ_b = (double) orb_occ;
    if (orb_spin == SPIN::Paired) occ_p = (double) orb_occ;
    
    double occ(0.0);
    if (dens_spin == DENSITY::Total) occ = occ_a + occ_b + occ_p;
    if (dens_spin == DENSITY::Alpha) occ = occ_a + 0.5*occ_p;
    if (dens_spin == DENSITY::Beta)  occ = occ_b + 0.5*occ_p;
    if (dens_spin == DENSITY::Spin)  occ = occ_a - occ_b;

    return occ;
}
    
} //namespace mrchem
