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
    densities.push_back(Density(false)); // d/dx rho total (meta-GGA, from orbital gradients)
    densities.push_back(Density(false)); // d/dy rho total
    densities.push_back(Density(false)); // d/dz rho total
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
            // grad rho = sum_i occ_i 2 phi_i grad phi_i from the SAME orbital gradients as tau (as in
            // MADNESS), so that tau_W = |grad rho|^2/(8 rho) and tau are consistent point by point:
            // for one orbital tau == tau_W exactly. With grad rho taken as the derivative of rho
            // instead, z = tau_W/tau scatters around its exact value wherever both are small
            // (e.g. one-orbital regions), and functionals that use tau_W/tau (TPSS, M06) amplify it.
            // Only for real orbitals; otherwise makepot() differentiates rho as before.
            bool real_orbs = true;
            for (auto &phi_i : *orbitals) real_orbs = real_orbs and (phi_i.iscomplex() == 0);
            Density *grad_rho[3] = {&densities[4], &densities[5], &densities[6]};
            if (tau.Ncomp() == 0) {
                auto derivOp = this->mrdft->functional().getTauDerivOp();
                if (derivOp == nullptr) MSG_ABORT("XCPotentialD1::setupDensities: derivative operator not set, cannot build tau");
                NablaOperator nabla(derivOp);
                nabla.setup(prec);

                // Same structure as density::compute(): each MPI rank sums the
                // contributions of its own orbitals, then the result is reduced.
                // Contributions are accumulated with the in-place CompFunction::add
                // (extends to union grid), never with mrcpp::add(out, ..., out, ...).
                // Products of orbital DERIVATIVES are built on the grid of the factors (maxIter = 0: no
                // adaptive refinement). A multiwavelet derivative is discontinuous at every node boundary,
                // so refining |d phi|^2 to prec subdivides along all of them: the XC grid (union of rho,
                // tau, grad rho) grew to ~30x the GGA size (H2 TPSS at 1e-6: 20000 vs 600 nodes, out of
                // memory).
                const int tau_maxiter = 0;
                Density tau_loc(false);
                std::vector<Density> grad_loc;
                for (int d = 0; d < 3; d++) grad_loc.push_back(Density(false));
                for (auto &phi_i : *orbitals) {
                    if (not mrcpp::mpi::my_func(phi_i)) continue;
                    double occ = static_cast<double>(phi_i.occ()); // total occupation, as for rho
                    if (std::abs(occ) < mrcpp::MachineZero) continue;
                    for (int d = 0; d < 3; d++) {
                        Orbital dphi = nabla[d](phi_i);
                        Density tau_d(false);
                        mrcpp::copy_grid(tau_d, dphi);
                        if (real_orbs) {
                            mrcpp::multiply(prec, tau_d, 1.0, dphi, dphi, tau_maxiter); // |d/dx_d phi_i|^2, no refinement
                        } else {
                            mrcpp::make_density(tau_d, dphi, prec); // |d/dx_d phi_i|^2 (real, also for complex phi)
                        }
                        tau_d.crop(prec);
                        tau_loc.add(0.5 * occ, tau_d);
                        if (real_orbs) { // d_d rho += 2 occ phi_i d_d phi_i
                            Density g_d(false);
                            mrcpp::copy_grid(g_d, dphi); // grid of the factors: phi's grid, widened by one node for a bandwidth-1 derivative
                            mrcpp::multiply(prec, g_d, 1.0, phi_i, dphi, tau_maxiter);
                            g_d.crop(prec);
                            grad_loc[d].add(2.0 * occ, g_d);
                        }
                    }
                }
                nabla.clear();
                density::allreduce_density(tau, tau_loc);
                if (real_orbs) {
                    for (int d = 0; d < 3; d++) density::allreduce_density(*grad_rho[d], grad_loc[d]);
                }
            }
            print_utils::qmfunction(3, "Compute tau", tau, timer);
            dens_vec.push_back(std::make_tuple(1.0, &tau.real()));
            if (grad_rho[0]->Ncomp() > 0) { // inputs 2, 3, 4: grad rho from orbital gradients (see makepot)
                for (int d = 0; d < 3; d++) dens_vec.push_back(std::make_tuple(1.0, &grad_rho[d]->real()));
            }
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
