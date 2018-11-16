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

#include "parallel.h"
#include "utils/math_utils.h"

#include "core.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace core {
/** @brief Helper struct to get the orbital ordering right
 *
 *  First index energy level (n)
 *  Second index angular momentum (l)
 */
int PT[29][2] = {
   /*s*/
   {1,0},                  /*p*/
   {2,0},                  {2,1},
   {3,0},            /*d*/ {3,1},
   {4,0},            {3,2},{4,1},
   {5,0},      /*f*/ {4,2},{5,1},
   {6,0},      {4,3},{5,2},{6,1},
   {7,0},/*g*/ {5,3},{6,2},{7,1},
   {8,0},{5,4},{6,3},{7,2},{8,1},
   {9,0},{6,4},{7,3},{8,2},{9,1}
};

} //namespace core
} //namespace initial_guess


/** @brief Produce an initial guess of orbitals
 *
 * @param prec: precision used in projection
 * @param mol: molecule
 * @param restricted: spin restriction
 * @param zeta: quality of hydrogen AO basis
 *
 * Sets up an AO basis of hydrogen functions with the given zeta quality
 * (SZ, DZ, TZ, QZ), computes and diagonalizes the core Hamiltonian matrix,
 * and fills the resulting orbitals by the Aufbau principle.
 *
 */
OrbitalVector initial_guess::core::setup(double prec,
                                         const Molecule &mol,
                                         bool restricted,
                                         int zeta) {
    int mult = mol.getMultiplicity();   //multiplicity
    int Ne = mol.getNElectrons();       //total electrons
    int Nd = Ne - (mult - 1);           //doubly occupied
    if (Nd%2 != 0) MSG_FATAL("Invalid multiplicity");

    // Make Fock operator contributions
    mrcpp::ABGVOperator<3> D(*MRA, 0.5, 0.5);
    KineticOperator T(D);
    NuclearOperator V(mol.getNuclei(), prec);

    // Project AO basis of hydrogen functions
    OrbitalVector Phi = initial_guess::core::project_ao(prec, mol.getNuclei(), SPIN::Paired, zeta);
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);

    // Compute Hamiltonian matrix
    Timer t_diag;
    Printer::printHeader(0, "Diagonalize Core-Hamiltonian matrix");
    Timer t1;
    T.setup(prec);
    V.setup(prec);
    ComplexMatrix t = T(Phi, Phi);
    ComplexMatrix v = V(Phi, Phi);
    ComplexMatrix F = S_m12.transpose()*(t + v)*S_m12;
    V.clear();
    T.clear();
    t1.stop();
    Printer::printDouble(0, "Compute Fock matrix", t1.getWallTime(), 5);

    // Diagonalize Hamiltonian matrix
    Timer t2;
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(F.cols());
    es.compute(F);
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexMatrix U = ei_vec.transpose() * S_m12;
    t2.stop();
    Printer::printDouble(0, "Diagonalize Fock matrix", t2.getWallTime(), 5);

    // Need to convert to QMFunctions for linear_combination
    QMFunctionVector funcs;
    for (auto &phi_i : Phi) funcs.push_back(phi_i);

    // Rotate orbitals and fill electrons by Aufbau
    Timer t3;
    OrbitalVector Psi;
    if (restricted) {
        if (mult != 1) MSG_FATAL("Restricted open-shell not available");

        int Np = Nd/2;                  //paired orbitals
        for (int i = 0; i < Np; i++) {
            ComplexVector v_i = U.row(i);
            Orbital psi_i(SPIN::Paired);
            qmfunction::linear_combination(psi_i, v_i, funcs, prec);
            Psi.push_back(psi_i);
        }
    } else {
        int Na = Nd/2 + (mult - 1);     //alpha orbitals
        int Nb = Nd/2;                  //beta orbitals

        OrbitalVector Psi_a;
        for (int i = 0; i < Na; i++) {
            ComplexVector v_i = U.row(i);
            Orbital psi_i(SPIN::Alpha);
            qmfunction::linear_combination(psi_i, v_i, funcs, prec);
            Psi_a.push_back(psi_i);
        }

        OrbitalVector Psi_b;
        for (int i = 0; i < Nb; i++) {
            ComplexVector v_i = U.row(i);
            Orbital psi_i(SPIN::Beta);
            qmfunction::linear_combination(psi_i, v_i, funcs, prec);
            Psi_b.push_back(psi_i);
        }

        Psi = orbital::adjoin(Psi_a, Psi_b);
    }
    orbital::free(Phi);
    t3.stop();
    Printer::printDouble(0, "Rotate orbitals", t3.getWallTime(), 5);

    t_diag.stop();
    Printer::printFooter(0, t_diag, 1);

    math_utils::print_matrix(0, es.eigenvalues(), "Eigenvalues", 10);

    return Psi;
}

/** @brief Project AO basis of hydrogen functions
 *
 * @param prec: precision used in projection
 * @param nucs: the nuclei of the molecule
 * @param zeta: quality of hydrogen AO basis
 *
 * Sets up an AO basis of hydrogen functions with the given zeta quality
 * (SZ, DZ, TZ, QZ), and projects it onto the MW basis. The basis at each
 * atomic center is always a complete shell: single zeta (SZ) means that
 * the current shell is completed, double zeta (DZ) means that also the
 * next shell will be completed, etc. E.i. the oxygen atom will get the
 * following AO basis:
 *
 * Oxygen AOs:
 * SZ: 1s2s2p                 (2s +  3p)
 * DZ: 1s2s2p3s3p             (3s +  6p)
 * TZ: 1s2s2p3s3p4s3d4p       (4s +  9p +  5d)
 * QZ: 1s2s2p3s3p4s3d4p5s4d5p (5s + 12p + 10d)
 *
 */
OrbitalVector initial_guess::core::project_ao(double prec,
                                              const Nuclei &nucs,
                                              int spin,
                                              int zeta) {
    Printer::printHeader(0, "Projecting Hydrogen AOs");
    println(0, "    N    Atom   Label                     SquareNorm");
    Printer::printSeparator(0, '-');

    Timer timer;
    OrbitalVector Phi;

    const char label[10] = "spdfg";

    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        int minAO = std::ceil(nuc.getElement().getZ()/2.0);
        double Z = nuc.getCharge();
        const mrcpp::Coord<3> &R = nuc.getCoord();

        int nAO = 0;
        int nShell = 0;
        int zetaReached = 0;
        bool minAOReached = false;
        while (true) {
            int n = initial_guess::core::PT[nShell][0];
            int l = initial_guess::core::PT[nShell][1];
            int M = 2*l + 1;

            if (minAOReached and l == 0) zetaReached++;
            if (zetaReached >= zeta) break;

            for (int m = 0; m < M; m++) {
                HydrogenFunction h_func(n, l, m, Z, R);

                Phi.push_back(spin);
                Phi.back().setRankID(Phi.size()%mpi::orb_size);
                Phi.back().alloc(NUMBER::Real);
                if (mpi::my_orb(Phi.back())) mrcpp::project(prec, Phi.back().real(), h_func);

                printout(0, std::setw(5)  << Phi.size());
                printout(0, std::setw(6)  << nuc.getElement().getSymbol() << i+1);
                printout(0, std::setw(6)  << n << label[l]);
                printout(0, std::setw(40) << Phi.back().squaredNorm());
                printout(0, std::endl);

                if (++nAO >= minAO) minAOReached = true;
            }
            nShell++;
        }
    }
    timer.stop();
    Printer::printFooter(0, timer, 2);
    return Phi;
}

} //namespace mrchem
