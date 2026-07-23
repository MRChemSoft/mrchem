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

        if (this->mrdft->functional().isMetaGGA()) { // Unperturbed kinetic energy density tau(r) = 1/2 sum_i |∇phi_i(r)|^2
            Timer timer;
            Density &tau = densities[3]; // tau total
            if (tau.Ncomp() == 0) {
                tau.alloc(1);
                mrcpp::copy_grid(tau.real(), grid);
            }
            tau.real().clear(); // start from zero

            // Use the derivative operator already stored in the functional
            auto derivOp = this->mrdft->functional().getDerivOp();
            if (!derivOp) {
                MSG_ABORT("XCPotentialD1::setupDensities: derivative operator not set, cannot build tau.\n");
            }

            // Temporary trees on the XC grid
            mrcpp::FunctionTree<3> tau_phi(grid.getMRA());
            mrcpp::copy_grid(tau_phi, grid);

            mrcpp::FunctionTree<3> d2phi_x(grid.getMRA());
            mrcpp::copy_grid(d2phi_x, grid);

            mrcpp::FunctionTree<3> d2phi_y(grid.getMRA());
            mrcpp::copy_grid(d2phi_y, grid);

            mrcpp::FunctionTree<3> d2phi_z(grid.getMRA());
            mrcpp::copy_grid(d2phi_z, grid);

            // Loop over orbitals: tau += 1/2 * tau_phi = 1/2 * |∇phi_i|^2
            for (int i = 0; i < static_cast<int>(orbitals->size()); ++i) {
                const mrcpp::CompFunction<3> &phi_i = (*orbitals)[i];
                mrcpp::FunctionTree<3> &phi_r =
                    const_cast<mrcpp::FunctionTree<3>&>(phi_i.real());

                tau_phi.clear();
                d2phi_x.clear();
                d2phi_y.clear();
                d2phi_z.clear();

                // Getting the gradient of orbitals
                mrcpp::FunctionTreeVector<3> grad_vec =
                    this->mrdft->functional().getLogGradient()
                        ? mrdft::xc_utils::log_gradient(*derivOp, phi_r)
                        : mrcpp::gradient(*derivOp, phi_r);

                mrcpp::FunctionTree<3> &dphi_x = mrcpp::get_func(grad_vec, 0);
                mrcpp::FunctionTree<3> &dphi_y = mrcpp::get_func(grad_vec, 1);
                mrcpp::FunctionTree<3> &dphi_z = mrcpp::get_func(grad_vec, 2);

                // mrcpp::FunctionTreeVector<3> hess_x_vec =
                //     this->mrdft->functional().getLogGradient()
                //         ? mrdft::xc_utils::log_gradient(*derivOp, dphi_x)
                //         : mrcpp::gradient(*derivOp, dphi_x);
                // mrcpp::FunctionTreeVector<3> hess_y_vec =
                //     this->mrdft->functional().getLogGradient()
                //         ? mrdft::xc_utils::log_gradient(*derivOp, dphi_y)
                //         : mrcpp::gradient(*derivOp, dphi_y);
                // mrcpp::FunctionTreeVector<3> hess_z_vec =
                //     this->mrdft->functional().getLogGradient()
                //         ? mrdft::xc_utils::log_gradient(*derivOp, dphi_z)
                //         : mrcpp::gradient(*derivOp, dphi_z);

                // mrcpp::FunctionTree<3> &d2phi_x = mrcpp::get_func(hess_x_vec, 0);
                // mrcpp::FunctionTree<3> &d2phi_y = mrcpp::get_func(hess_y_vec, 1);
                // mrcpp::FunctionTree<3> &d2phi_z = mrcpp::get_func(hess_z_vec, 2);


                // mrcpp::add(prec, tau_phi, 1.0, tau_phi, 1.0, d2phi_x);  // tau_phi += d2phi_x
                // mrcpp::add(prec, tau_phi, 1.0, tau_phi, 1.0, d2phi_y);  // tau_phi += d2phi_y
                // mrcpp::add(prec, tau_phi, 1.0, tau_phi, 1.0, d2phi_z);  // tau_phi += d2phi_z


                // Constructing tau_phi as sum_i=xyz dphi_i^2
                mrcpp::power(prec, d2phi_x, dphi_x, 2.0);              // d2phi_x = (dphi_x)^2
                mrcpp::power(prec, d2phi_y, dphi_y, 2.0);              // d2phi_y = (dphi_y)^2
                mrcpp::power(prec, d2phi_z, dphi_z, 2.0);              // d2phi_z = (dphi_z)^2

                mrcpp::add(prec, tau_phi, 1.0, tau_phi, 1.0, d2phi_x);  // d2phi_ += d2phi_x
                mrcpp::add(prec, tau_phi, 1.0, tau_phi, 1.0, d2phi_y);  // d2phi_ += d2phi_y
                mrcpp::add(prec, tau_phi, 1.0, tau_phi, 1.0, d2phi_z);  // d2phi_ += d2phi_z

                // tau needs the same orbital occupancy weight as the density.
                // For paired orbitals this restores the usual factor of 2.
                double occ = phi_i.occ();
                mrcpp::add(prec, tau.real(), 1.0, tau.real(), 0.5 * occ, tau_phi);
            }

            print_utils::qmfunction(3, "Compute tau", tau, timer);
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
