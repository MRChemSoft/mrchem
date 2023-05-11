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

#include "CoulombPotential.h"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

/** @brief constructor
 *
 * @param P: MW Poisson operator
 * @param Phi: orbitals defining the operator
 *
 * Density is not spin-separated since this operator requires only total density.
 * This operator will always point to the same OrbitalVector, but the orbitals within
 * the vector can change throughout the calculation. The density and (*this)
 * QMPotential is uninitialized at this point and will be computed at setup.
 */

CoulombPotential::CoulombPotential(PoissonOperator_p P, OrbitalVector_p Phi, bool mpi_share)
        : QMPotential(1, mpi_share)
        , density(false)
        , orbitals(Phi)
        , poisson(P) {}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application of the Poisson
 * operator to the density. If the density is not available it is computed
 * from the current orbitals (assuming that the orbitals are available).
 * For first-order perturbations the first order density and the Hessian will be
 * computed. In order to make the Hessian available to CoulombOperator, it is stored in the
 * potential function instead of the zeroth-order potential.
 *
 */
void CoulombPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    Timer timer;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(3, "Building Coulomb operator");
    mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
    mrcpp::print::separator(3, '-');
    if (hasDensity()) {
        setupGlobalPotential(prec);
    } else if (mrcpp::mpi::numerically_exact) {
        setupGlobalDensity(prec);
        setupGlobalPotential(prec);
    } else {
        // Keep each local contribution a bit
        // more precise than strictly necessary
        setupLocalDensity(0.1 * prec);
        mrcpp::ComplexFunction V = setupLocalPotential(0.1 * prec);
        allreducePotential(0.1 * prec, V);
    }
    if (plevel == 2) print_utils::qmfunction(2, "Coulomb operator", *this, timer);
    mrcpp::print::footer(3, timer, 2);
}

/** @brief clear operator after application
 *
 * This will clear the operator and bring it back to the state after construction.
 * The operator can now be reused after another setup.
 */
void CoulombPotential::clear() {
    mrcpp::ComplexFunction::free(NUMBER::Total); // delete FunctionTree pointers
    this->density.free(NUMBER::Total);           // delete FunctionTree pointers
    clearApplyPrec();                            // apply_prec = -1
}

/** @brief compute Coulomb potential
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application o the Poisson operator
 * to the precomputed electron density.
 */
void CoulombPotential::setupGlobalPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    mrcpp::ComplexFunction &V = *this;
    mrcpp::ComplexFunction &rho = this->density;

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    // Adjust precision by system size
    double abs_prec = prec / rho.norm();
    bool need_to_apply = not(V.isShared()) or mrcpp::mpi::share_master();

    Timer timer;
    V.alloc(NUMBER::Real);
    if (need_to_apply) mrcpp::apply(abs_prec, V.real(), P, rho.real());
    mrcpp::mpi::share_function(V, 0, 22445, mrcpp::mpi::comm_share);
    print_utils::qmfunction(3, "Compute global potential", V, timer);
}

/** @brief compute Coulomb potential
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application o the Poisson operator
 * to the precomputed electron density.
 */
mrcpp::ComplexFunction CoulombPotential::setupLocalPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    OrbitalVector &Phi = *this->orbitals;
    mrcpp::ComplexFunction &rho = this->density;

    // Adjust precision by system size
    double abs_prec = prec / orbital::get_electron_number(Phi);

    Timer timer;
    mrcpp::ComplexFunction V(false);
    V.alloc(NUMBER::Real);
    mrcpp::apply(abs_prec, V.real(), P, rho.real());
    print_utils::qmfunction(3, "Compute local potential", V, timer);

    return V;
}

void CoulombPotential::allreducePotential(double prec, mrcpp::ComplexFunction &V_loc) {
    Timer t_com;

    mrcpp::ComplexFunction &V_tot = *this;
    OrbitalVector &Phi = *this->orbitals;

    double abs_prec = prec / orbital::get_electron_number(Phi);

    // Add up local contributions into the grand master
    mrcpp::mpi::reduce_function(abs_prec, V_loc, mrcpp::mpi::comm_wrk);

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
    print_utils::qmfunction(3, "Allreduce potential", V_tot, t_com);
}

} // namespace mrchem
