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

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "XCPotentialD1.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include "mrdft/xc_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

XCPotentialD1::XCPotentialD1(std::unique_ptr<mrdft::MRDFT> &F, std::shared_ptr<OrbitalVector> Phi, bool mpi_shared)
        : XCPotential(F, Phi, mpi_shared) {
    densities.push_back(Density(false)); // rho_0 total
    densities.push_back(Density(false)); // rho_0 alpha
    densities.push_back(Density(false)); // rho_0 beta
    densities.push_back(Density(false)); // tau total
    // densities.push_back(Density(false)); // tau alpha
    // densities.push_back(Density(false)); // tau beta
}

mrcpp::FunctionTreeVector<3> XCPotentialD1::setupDensities(double prec, mrcpp::FunctionTree<3> &grid) {
    mrcpp::FunctionTreeVector<3> dens_vec;
    if (not this->mrdft->functional().isSpin()) {

        { // Unperturbed total density
            Timer timer;
            Density &rho = getDensity(DensityType::Total, 0);
            if (rho.Ncomp() == 0) {
                rho.alloc(1);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DensityType::Total);
            }
            print_utils::qmfunction(3, "Compute rho", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }

        if (this->mrdft->functional().isMetaGGA()) { // Unperturbed kinetic energy density tau(r) = 1/2 sum_i occ_i |∇phi_i(r)|^2
            Timer timer;
            // NB: tau is stored in densities[3]; do not use getDensity(Total, 1) here,
            // that index is reserved for the perturbed density in the response code.
            Density &tau = densities[3];
            if (tau.Ncomp() == 0) {
                auto derivOp = this->mrdft->functional().getDerivOp();
                if (derivOp == nullptr) MSG_ABORT("XCPotentialD1::setupDensities: derivative operator not set, cannot build tau");
                NablaOperator nabla(derivOp);
                nabla.setup(prec);

                // Same structure as density::compute(): each MPI rank sums the
                // contributions of its own orbitals, then the result is reduced.
                // Contributions are accumulated with the in-place CompFunction::add
                // (extends to union grid), never with mrcpp::add(out, ..., out, ...).
                Density tau_loc(false);
                for (auto &phi_i : *orbitals) {
                    if (not mrcpp::mpi::my_func(phi_i)) continue;
                    double occ = static_cast<double>(phi_i.occ()); // total occupation, as for rho
                    if (std::abs(occ) < mrcpp::MachineZero) continue;
                    for (int d = 0; d < 3; d++) {
                        Orbital dphi = nabla[d](phi_i);
                        Density tau_d(false);
                        mrcpp::copy_grid(tau_d, dphi);
                        mrcpp::make_density(tau_d, dphi, prec); // |d/dx_d phi_i|^2 (real, also for complex phi)
                        tau_d.crop(prec);
                        tau_loc.add(0.5 * occ, tau_d);
                    }
                }
                nabla.clear();
                density::allreduce_density(tau, tau_loc);
            }
            print_utils::qmfunction(3, "Compute tau", tau, timer);
            // Debug check: for real orbitals int tau = sum_i occ_i <phi_i|T|phi_i> = E_kin
            mrcpp::print::value(1, "Integral of tau (=E_kin)", tau.integrate().real(), "(au)", 10);
            dens_vec.push_back(std::make_tuple(1.0, &tau.real()));
        }
    } else {
        { // Unperturbed alpha density
            Timer timer;
            Density &rho = getDensity(DensityType::Alpha, 0);
            if (rho.Ncomp() == 0) {
                rho.alloc(1);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DensityType::Alpha);
            }
            print_utils::qmfunction(3, "Compute rho (alpha)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
        { // Unperturbed beta density
            Timer timer;
            Density &rho = getDensity(DensityType::Beta, 0);
            if (rho.Ncomp() == 0) {
                rho.alloc(1);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DensityType::Beta);
            }
            print_utils::qmfunction(3, "Compute rho (beta)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }

        if (this->mrdft->functional().isMetaGGA()) {
            MSG_ABORT("Spin-polarized meta-GGA not supported in XCPotentialD1::setupDensities.");
        }
    }
    return dens_vec;
}

} // namespace mrchem
