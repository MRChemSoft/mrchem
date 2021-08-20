/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/MWOperators"
#include "MRCPP/Plotter"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "FarFieldPotential.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

#include <sstream>
using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

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

FarFieldPotential::FarFieldPotential(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, const Nuclei &nucs, double prec, bool mpi_share)
        : QMPotential(1, mpi_share)
        , exp_prec(prec)
        , nuclei(nucs)
        , density(false)
        , orbitals(Phi)
        , poisson(P) {}

void FarFieldPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    if (hasDensity()) {
        setupPotential(prec);
    } else {
        setupDensity(prec);
        setupPotential(prec);
    }
}

void FarFieldPotential::setupPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    QMFunction &V = *this;
    QMFunction &rho = this->density;

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    // Adjust precision by system size
    double abs_prec = prec / rho.norm();
    // double abs_prec = prec;
    bool need_to_apply = not(V.isShared()) or mpi::share_master();
    Timer timer;
    V.alloc(NUMBER::Real);
    if (need_to_apply) {
        if ((*MRA).getOperatorScale() == 0) {
            mrcpp::apply_far_field(abs_prec, V.real(), P, rho.real());
        } else if ((*MRA).getOperatorScale() < 0) {
            // Old implementation
            mrcpp::FunctionTree<3> near_field(*MRA);
            mrcpp::FunctionTree<3> full_field(*MRA);
            mrcpp::apply_near_field(abs_prec, near_field, P, rho.real());
            mrcpp::apply(abs_prec, full_field, P, rho.real());
            mrcpp::build_grid(V.real(), near_field);
            mrcpp::build_grid(V.real(), full_field);
            mrcpp::add(-1, V.real(), 1.0, full_field, -1.0, near_field);
        } else {
            MSG_ABORT("Operator scale cannot be positive");
        }
    }

    mpi::share_function(V, 0, 22445, mpi::comm_share);
    print_utils::qmfunction(2, "FarField potential", V, timer);
}

void FarFieldPotential::setupDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;

    Density &rho = this->density;

    Density rho_el(false);
    density::compute(prec, rho_el, Phi, DensityType::Total);
    auto rho_nuc = chemistry::compute_nuclear_density(prec, this->nuclei, this->exp_prec);
    qmfunction::add(rho, 1.0, rho_el, -1.0, rho_nuc, -1.0);
    if (std::norm(rho.integrate()) > this->apply_prec) MSG_ABORT("Non-zero net charge on unit cell");
}

void FarFieldPotential::clear() {
    QMFunction::free(NUMBER::Total);   // delete FunctionTree pointers
    this->density.free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();                  // apply_prec = -1
}

} // namespace mrchem
