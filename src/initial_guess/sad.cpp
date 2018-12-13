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

#include <Eigen/Eigenvalues>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "core.h"
#include "gto.h"
#include "parallel.h"
#include "sad.h"
#include "utils/math_utils.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace sad {

ComplexMatrix diagonalize_fock(RankZeroTensorOperator &F, OrbitalVector &Phi, int spin);
OrbitalVector rotate_orbitals(double prec, ComplexMatrix &U, OrbitalVector &Phi, int N, int spin);
void project_atomic_densities(double prec, const Molecule &mol, mrcpp::FunctionTreeVector<3> &rho_atomic);

} //namespace sad
} //namespace initial_guess

OrbitalVector initial_guess::sad::setup(double prec, const Molecule &mol, bool restricted, int zeta) {
    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity(); //multiplicity
    int Ne = mol.getNElectrons();     //total electrons
    int Nd = Ne - (mult - 1);         //doubly occupied electrons
    if (Nd % 2 != 0) MSG_FATAL("Invalid multiplicity");
    int Na = Nd / 2 + (mult - 1); //alpha orbitals
    int Nb = Nd / 2;              //beta orbitals

    // Make Fock operator contributions
    mrcpp::PoissonOperator P(*MRA, prec);
    mrcpp::ABGVOperator<3> D(*MRA, 0.5, 0.5);
    mrdft::XCFunctional xcfun(*MRA, not restricted);
    xcfun.setFunctional("SLATERX");
    xcfun.setFunctional("VWN5C");
    xcfun.evalSetup(1);
    KineticOperator T(D);
    NuclearOperator V(mol.getNuclei(), prec);
    CoulombOperator J(&P);
    XCOperator XC(&xcfun);
    RankZeroTensorOperator F = T + V + J + XC;

    // Compute Coulomb density
    Density &rho_j = J.getDensity();
    rho_j.function().alloc(NUMBER::Real);

    // MPI grand master computes SAD density
    if (mpi::grand_master()) {
        // Compute atomic densities
        mrcpp::FunctionTreeVector<3> rho_atomic;
        initial_guess::sad::project_atomic_densities(prec, mol, rho_atomic);

        // Add atomic densities
        mrcpp::add(prec, rho_j.function().real(), rho_atomic);
        mrcpp::clear(rho_atomic, true);
    }
    // MPI grand master distributes the full density
    mpi::broadcast_density(rho_j, mpi::comm_orb);

    // Compute XC density
    if (restricted) {
        mrcpp::FunctionTree<3> &rho_xc = XC.getDensity(DENSITY::Total);
        mrcpp::copy_grid(rho_xc, rho_j.function().real());
        mrcpp::copy_func(rho_xc, rho_j.function().real());
    } else {
        mrcpp::FunctionTree<3> &rho_a = XC.getDensity(DENSITY::Alpha);
        mrcpp::FunctionTree<3> &rho_b = XC.getDensity(DENSITY::Beta);
        mrcpp::add(prec, rho_a, 1.0, rho_j.function().real(), -1.0 * Nb / Ne, rho_j.function().real());
        mrcpp::add(prec, rho_b, 1.0, rho_j.function().real(), -1.0 * Na / Ne, rho_j.function().real());

        // Extend to union grid
        int nNodes = 1;
        while (nNodes > 0) {
            int nAlpha = mrcpp::refine_grid(rho_a, rho_b);
            int nBeta = mrcpp::refine_grid(rho_b, rho_a);
            nNodes = nAlpha + nBeta;
        }
    }

    // Project AO basis of hydrogen functions
    OrbitalVector Phi = initial_guess::core::project_ao(prec, mol.getNuclei(), SPIN::Paired, zeta);

    Timer t_fock;
    Printer::printHeader(0, "Setting up Fock operator");
    F.setup(prec);
    t_fock.stop();
    Printer::printFooter(0, t_fock, 2);

    // Compute Fock matrix
    Timer t_diag;
    Printer::printHeader(0, "Diagonalize Fock matrix");
    OrbitalVector Psi;
    if (restricted) {
        if (mult != 1) MSG_FATAL("Restricted open-shell not available");
        int Np = Nd / 2; //paired orbitals
        ComplexMatrix U = initial_guess::sad::diagonalize_fock(F, Phi, SPIN::Paired);
        Psi = initial_guess::sad::rotate_orbitals(prec, U, Phi, Np, SPIN::Paired);
    } else {
        int Na = Nd / 2 + (mult - 1); //alpha orbitals
        int Nb = Nd / 2;              //beta orbitals

        ComplexMatrix U_a = initial_guess::sad::diagonalize_fock(F, Phi, SPIN::Alpha);
        OrbitalVector Psi_a = initial_guess::sad::rotate_orbitals(prec, U_a, Phi, Na, SPIN::Alpha);

        ComplexMatrix U_b = initial_guess::sad::diagonalize_fock(F, Phi, SPIN::Beta);
        OrbitalVector Psi_b = initial_guess::sad::rotate_orbitals(prec, U_b, Phi, Nb, SPIN::Beta);

        Psi = orbital::adjoin(Psi_a, Psi_b);
    }
    F.clear();
    t_diag.stop();
    Printer::printFooter(0, t_diag, 2);

    return Psi;
}

ComplexMatrix initial_guess::sad::diagonalize_fock(RankZeroTensorOperator &F, OrbitalVector &Phi, int spin) {
    Timer t1;
    for (int i = 0; i < Phi.size(); i++) Phi[i].setSpin(spin);
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    ComplexMatrix f_tilde = F(Phi, Phi);
    ComplexMatrix f = S_m12.adjoint() * f_tilde * S_m12;
    t1.stop();
    Printer::printDouble(0, "Compute Fock matrix", t1.getWallTime(), 5);

    Timer t2;
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(f.cols());
    es.compute(f);
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexMatrix U = ei_vec.transpose() * S_m12;
    t2.stop();
    Printer::printDouble(0, "Diagonalize Fock matrix", t2.getWallTime(), 5);
    return U;
}

OrbitalVector initial_guess::sad::rotate_orbitals(double prec, ComplexMatrix &U, OrbitalVector &Phi, int N, int spin) {
    Timer t;
    OrbitalVector Psi;
    for (int i = 0; i < N; i++) Psi.push_back(Orbital(spin));
    mpi::distribute(Psi);

    OrbitalIterator iter(Phi);
    while (iter.next()) {
        for (int i = 0; i < Psi.size(); i++) {
            if (not mpi::my_orb(Psi[i])) continue;
            QMFunctionVector func_vec;
            ComplexVector coef_vec(iter.get_size());
            for (int j = 0; j < iter.get_size(); j++) {
                int idx_j = iter.idx(j);
                Orbital &recv_j = iter.orbital(j);
                coef_vec[j] = U(i, idx_j);
                func_vec.push_back(recv_j);
            }
            Orbital tmp_i = Psi[i].paramCopy();
            qmfunction::linear_combination(tmp_i, coef_vec, func_vec, prec);
            Psi[i].add(1.0, tmp_i); // In place addition
            Psi[i].crop(prec);
        }
    }
    t.stop();
    Printer::printDouble(0, "Rotate orbitals", t.getWallTime(), 5);
    return Psi;
}

void initial_guess::sad::project_atomic_densities(double prec,
                                                  const Molecule &mol,
                                                  mrcpp::FunctionTreeVector<3> &rho_atomic) {
    Timer timer;
    Printer::printHeader(0, "Projecting Gaussian-type density");
    println(0, " Nr  Element                                 Rho_i");
    Printer::printSeparator(0, '-');

    std::string sad_path = SAD_BASIS_DIR;

    int oldprec = Printer::setPrecision(15);
    const Nuclei &nucs = mol.getNuclei();
    for (int k = 0; k < nucs.size(); k++) {
        const std::string &sym = nucs[k].getElement().getSymbol();

        std::stringstream bas;
        std::stringstream dens;
        bas << sad_path << sym << ".bas";
        dens << sad_path << sym << ".dens";

        mrcpp::FunctionTree<3> *rho = initial_guess::gto::project_density(prec, nucs[k], bas.str(), dens.str());
        printout(0, std::setw(3) << k);
        printout(0, std::setw(7) << sym);
        printout(0, std::setw(49) << rho->integrate() << "\n");

        rho_atomic.push_back(std::make_tuple(1.0, rho));
    }
    Printer::setPrecision(oldprec);

    timer.stop();
    Printer::printFooter(0, timer, 2);
}

} //namespace mrchem
