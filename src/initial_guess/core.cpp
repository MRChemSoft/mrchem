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

#include <MRCPP/MWOperators>
#include <MRCPP/Parallel>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "core.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "chemistry/Nucleus.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/MomentumOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/qmoperator_utils.h"

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
// clang-format off
int PT[29][2] = {
     /*s*/
    {1, 0},                      /*p*/
    {2, 0},                     {2, 1},
    {3, 0},               /*d*/ {3, 1},
    {4, 0},              {3, 2},{4, 1},
    {5, 0},        /*f*/ {4, 2},{5, 1},
    {6, 0},       {4, 3},{5, 2},{6, 1},
    {7, 0}, /*g*/ {5, 3},{6, 2},{7, 1},
    {8, 0},{5, 4},{6, 3},{7, 2},{8, 1},
    {9, 0},{6, 4},{7, 3},{8, 2},{9, 1}
};
// clang-format on

} // namespace core
} // namespace initial_guess

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
bool initial_guess::core::setup(OrbitalVector &Phi, double prec, const Nuclei &nucs, int zeta) {
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Compute initial orbitals");
    print_utils::text(0, "Method      ", "Diagonalize Core Hamiltonian matrix");
    print_utils::text(0, "Precision   ", print_utils::dbl_to_str(prec, 5, true));
    print_utils::text(0, "Restricted  ", (orbital::size_singly(Phi)) ? "False" : "True");
    print_utils::text(0, "AO basis    ", "Hydrogenic orbitals");
    print_utils::text(0, "Zeta quality", std::to_string(zeta));
    mrcpp::print::separator(0, '~', 2);

    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    plevel = 1;
    if (plevel == 1) mrcpp::print::header(1, "Core-Hamiltonian Initial Guess");

    // Make Fock operator contributions
    t_lap.start();
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
    MomentumOperator p(D_p);
    NuclearOperator V(nucs, prec);
    if (plevel == 1) mrcpp::print::time(1, "Projecting nuclear potential", t_lap);

    // Project AO basis of hydrogen functions
    t_lap.start();
    OrbitalVector Psi;
    initial_guess::core::project_ao(Psi, prec, nucs, zeta);
    if (plevel == 1) mrcpp::print::time(1, "Projecting Hydrogen AOs", t_lap);

    p.setup(prec);
    V.setup(prec);

    // Compute Hamiltonian matrix
    t_lap.start();
    mrcpp::print::header(2, "Diagonalize Hamiltonian matrix");
    ComplexMatrix U = initial_guess::core::diagonalize(Psi, p, V);

    // Rotate orbitals and fill electrons by Aufbau
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);

    initial_guess::core::rotate_orbitals(Phi, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_a, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_b, prec, U, Psi);
    Phi = orbital::adjoin(Phi, Phi_a);
    Phi = orbital::adjoin(Phi, Phi_b);

    V.clear();
    p.clear();

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
    return true;
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
void initial_guess::core::project_ao(OrbitalVector &Phi, double prec, const Nuclei &nucs, int zeta) {
    Timer t_tot;
    auto w0 = Printer::getWidth() - 2;
    auto w1 = 5;
    auto w2 = 7;
    auto w3 = w0 * 2 / 9;
    auto w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w4) << "Atom";
    o_head << std::setw(w2) << "Label";
    o_head << std::setw(w3 + 1) << "Nodes";
    o_head << std::setw(w3) << "Size";
    o_head << std::setw(w3) << "Time";

    mrcpp::print::header(2, "Projecting Hydrogen AOs");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    const char label[10] = "spdfg";

    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        int minAO = std::ceil(nuc.getElement().getZ() / 2.0);
        double Z = nuc.getCharge();
        const mrcpp::Coord<3> &R = nuc.getCoord();

        int nAO = 0;
        int nShell = 0;
        int zetaReached = 0;
        bool minAOReached = false;
        while (true) {
            int n = initial_guess::core::PT[nShell][0];
            int l = initial_guess::core::PT[nShell][1];
            int M = 2 * l + 1;

            if (minAOReached and l == 0) zetaReached++;
            if (zetaReached >= zeta) break;

            for (int m = 0; m < M; m++) {
                Timer t_i;
                HydrogenFunction h_func(n, l, m, Z, R);
                Phi.push_back(Orbital(SPIN::Paired));
                Phi.back().setRank(Phi.size() - 1);
                if (mrcpp::mpi::my_orb(Phi.back())) {
                    mrcpp::cplxfunc::project(Phi.back(), h_func, NUMBER::Real, prec);
                    if (std::abs(Phi.back().norm() - 1.0) > 0.01) MSG_WARN("AO not normalized!");
                }

                std::stringstream o_txt;
                o_txt << std::setw(w1 - 1) << Phi.size() - 1;
                o_txt << std::setw(w4) << nuc.getElement().getSymbol();
                o_txt << std::setw(w2 - 1) << n << label[l];
                print_utils::qmfunction(2, o_txt.str(), Phi.back(), t_i);

                if (++nAO >= minAO) minAOReached = true;
            }
            nShell++;
        }
    }
    mrcpp::print::footer(2, t_tot, 2);
}

void initial_guess::core::rotate_orbitals(OrbitalVector &Psi, double prec, ComplexMatrix &U, OrbitalVector &Phi) {
    if (Psi.size() == 0) return;
    Timer t_tot;
    mrcpp::mpifuncvec::rotate(Phi, U, Psi, prec);
    mrcpp::print::time(1, "Rotating orbitals", t_tot);
}

ComplexMatrix initial_guess::core::diagonalize(OrbitalVector &Phi, MomentumOperator &p, RankZeroOperator &V) {
    Timer t1;
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    mrcpp::print::separator(2, '-');
    ComplexMatrix t_tilde = qmoperator::calc_kinetic_matrix(p, Phi, Phi);
    ComplexMatrix v_tilde = V(Phi, Phi);
    ComplexMatrix f_tilde = t_tilde + v_tilde;
    ComplexMatrix f = S_m12.adjoint() * f_tilde * S_m12;
    mrcpp::print::separator(2, '-');
    mrcpp::print::time(1, "Computing Fock matrix", t1);

    Timer t2;
    DoubleVector eig;
    ComplexMatrix U = math_utils::diagonalize_hermitian_matrix(f, eig);
    mrcpp::print::time(1, "Diagonalizing Fock matrix", t2);

    return S_m12 * U;
}

} // namespace mrchem
