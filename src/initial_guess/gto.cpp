/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <MRCPP/Gaussians>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "gto.h"
#include "parallel.h"

#include "utils/gto_utils/Intgrl.h"
#include "utils/gto_utils/OrbitalExp.h"
#include "utils/math_utils.h"
#include "utils/print_utils.h"

#include "chemistry/Nucleus.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::GaussExp;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Produce an initial guess of orbitals
 *
 * @param Phi: vector or MW orbitals
 * @param prec: precision used in projection
 * @param bas_file: basis set file (LSDalton format)
 * @param mop_file: file with paired MO coefficients
 * @param moa_file: file with alpha MO coefficients
 * @param mob_file: file with beta MO coefficients
 *
 * Sets up a precomputed MO basis from an unrestricted LSDalton calculation.
 * Requires the LSDalton basis file and the corresponding MO matrices (not in
 * any official format!). The MO files should start with one entry giving the
 * number of AOs, followed by the columns of the MO matrix concatenated into
 * a single column.
 *
 * Projects only the occupied orbitals of each spin.
 *
 */
bool initial_guess::gto::setup(OrbitalVector &Phi,
                               double prec,
                               const std::string &bas_file,
                               const std::string &mop_file,
                               const std::string &moa_file,
                               const std::string &mob_file) {
    if (Phi.size() == 0) return false;

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation   ", "Read orbitals from file (GTO)");
    print_utils::text(0, "Precision     ", print_utils::dbl_to_str(prec, 5, true));
    if (orbital::size_singly(Phi)) {
        print_utils::text(0, "Restricted    ", "False");
        print_utils::text(0, "MO alpha file ", moa_file);
        print_utils::text(0, "MO beta file  ", mob_file);
    } else {
        print_utils::text(0, "Restricted    ", "True");
        print_utils::text(0, "MO file ", mop_file);
    }
    mrcpp::print::separator(0, '~', 2);

    // Separate alpha/beta from paired orbitals
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);

    // Project paired, alpha and beta separately
    initial_guess::gto::project_mo(Phi, prec, bas_file, mop_file);
    initial_guess::gto::project_mo(Phi_a, prec, bas_file, moa_file);
    initial_guess::gto::project_mo(Phi_b, prec, bas_file, mob_file);

    // Collect orbitals into one vector
    for (auto &phi_a : Phi_a) Phi.push_back(phi_a);
    for (auto &phi_b : Phi_b) Phi.push_back(phi_b);

    return true;
}

/** @brief Project the N first GTO expansions of the MO basis
 *
 * @param Phi: vector or MW orbitals
 * @param prec Precision used in projection
 * @param bas_file String containing basis set file
 * @param mo_file String containing MO matrix file
 *
 * Projects the N first rows of the MO matrix from GTO orbitals into
 * corresponding MW orbitals.
 *
 */
void initial_guess::gto::project_mo(OrbitalVector &Phi,
                                    double prec,
                                    const std::string &bas_file,
                                    const std::string &mo_file) {
    if (Phi.size() == 0) return;

    Timer t_tot;
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 2;
    auto w1 = 5;
    auto w2 = w0 * 2 / 9;
    auto w3 = w0 - w1 - 3 * w2;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w3) << "Norm";
    o_head << std::setw(w2 + 1) << "Nodes";
    o_head << std::setw(w2) << "Size";
    o_head << std::setw(w2) << "Time";

    mrcpp::print::header(1, "GTO Initial Guess");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    // Setup AO basis
    Timer t1;
    gto_utils::Intgrl intgrl(bas_file);
    gto_utils::OrbitalExp gto_exp(intgrl);
    t1.stop();

    // Read MO file (transpose)
    Timer t2;
    DoubleMatrix MO = math_utils::read_matrix_file(mo_file);
    if (MO.cols() < Phi.size()) MSG_ABORT("Size mismatch");
    t2.stop();

    Timer t3;
    for (int i = 0; i < Phi.size(); i++) {
        Timer t_i;
        if (mpi::my_orb(Phi[i])) {
            GaussExp<3> mo_i = gto_exp.getMO(i, MO.transpose());
            Phi[i].alloc(NUMBER::Real);
            mrcpp::project(prec, Phi[i].real(), mo_i);
        }
        std::stringstream o_txt;
        o_txt << std::setw(w1 - 1) << i;
        o_txt << std::setw(w3) << print_utils::dbl_to_str(Phi[i].norm(), pprec, true);
        print_utils::qmfunction(2, o_txt.str(), Phi[i], t_i);
    }
    mpi::barrier(mpi::comm_orb);
    mrcpp::print::separator(2, '-');
    mrcpp::print::time(1, "Reading AO basis", t1);
    mrcpp::print::time(1, "Reading MO matrix", t2);
    mrcpp::print::time(1, "Projecting GTO MOs", t3);
    mrcpp::print::footer(1, t_tot, 2);
}

/** @brief Project the N first GTO expansions of the AO basis
 *
 * @param Phi: vector or MW orbitals
 * @param prec Precision used in projection
 * @param bas_file String containing basis set file
 *
 * Projects the N first Gaussian-type AOs into corresponding MW orbitals.
 *
 */
void initial_guess::gto::project_ao(OrbitalVector &Phi, double prec, const std::string &bas_file) {
    if (Phi.size() == 0) return;
    mrcpp::print::header(0, "Setting up Gaussian-type AOs");
    println(1, "    n  Spin  Occ                           SquareNorm");
    mrcpp::print::separator(0, '-');
    Timer timer;

    // Setup AO basis
    gto_utils::Intgrl intgrl(bas_file);
    gto_utils::OrbitalExp gto_exp(intgrl);
    if (gto_exp.size() < Phi.size()) MSG_ABORT("Size mismatch");

    for (int i = 0; i < Phi.size(); i++) {
        GaussExp<3> ao_i = gto_exp.getAO(i);
        qmfunction::project(Phi[i], ao_i, NUMBER::Real, prec);
        printout(0, std::setw(5) << i);
        printout(0, std::setw(5) << Phi[i].printSpin());
        printout(0, std::setw(5) << Phi[i].occ());
        printout(0, std::setw(44) << Phi[i].norm() << std::endl);
        Phi.push_back(Phi[i]);
    }
    mpi::barrier(mpi::comm_orb);
    timer.stop();
    mrcpp::print::footer(0, timer, 2);
}

Density initial_guess::gto::project_density(double prec,
                                            const Nucleus &nuc,
                                            const std::string &bas_file,
                                            const std::string &dens_file) {
    // Setup AO basis
    gto_utils::Intgrl intgrl(bas_file);
    intgrl.getNucleus(0).setCoord(nuc.getCoord());
    gto_utils::OrbitalExp gto_exp(intgrl);

    // Read density matrix file
    DoubleMatrix D = math_utils::read_matrix_file(dens_file);
    GaussExp<3> dens_exp = gto_exp.getDens(D);

    Density rho(false);
    density::compute(prec, rho, dens_exp);
    return rho;
}

} // namespace mrchem

// void OrbitalVector::readVirtuals(const string &bf, const string &mo, int n_occ) {
//    Timer timer;
//    int oldPrec = Printer::setPrecision(15);
//    printout(0, "\n\n=============== Setting up virtual orbitals ");
//    printout(0, "================\n\n");

//    OrbitalExp *moExp = readOrbitalExpansion(bf, mo);
//    for (int a = n_occ; a < moExp->size(); a++) {
//    GaussExp<3> &gtOrb = moExp->getOrbital(a);
//        Orbital *orb_a = new Orbital(2, Orbital::Paired);
//    orb_a->projectFunction(gtOrb);
//        printout(0, "Orbital " << setw(3) << a);
//        println(0, " squareNorm: " << setw(36) << orb_a->getSquareNorm());
//        this->orbitals.push_back(orb_a);
//    }
//    delete moExp;
//    Printer::setPrecision(5);
//    printout(0, "\n================ Elapsed time: ");
//    println(0, timer.elapsed() << " =================\n");
//    Printer::setPrecision(oldPrec);
//}
