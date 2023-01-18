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

#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <MRCPP/utils/details.h>

#include "SCFSolver.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Set convergence thresholds
 *
 * @param orb: threshold for orbitals
 * @param prop: threshold for property
 */
void SCFSolver::setThreshold(double orb, double prop) {
    this->orbThrs = orb;
    this->propThrs = prop;
}

/** @brief Set dynamic precision parameters
 *
 * @param init: initial precision
 * @param final: final precision
 *
 * The precision will increase dynamically during the SCF optimization, starting
 * from "init", ending at "final".
 */
void SCFSolver::setOrbitalPrec(double init, double final) {
    this->orbPrec[0] = init;
    this->orbPrec[1] = init;
    this->orbPrec[2] = final;
}

/** @brief Reset accumulated data */
void SCFSolver::reset() {
    this->error.clear();
    this->property.clear();
    this->orbPrec[0] = this->orbPrec[1];
}

/** @brief Adjust dynamic precision
 *
 * @param error: error in current SCF iteration
 *
 * The precision will increase dynamically following the current residual in
 * the SCF equations, and at least by 25%. Sets the internal precision parameter
 * and returns the current prec.
 */
double SCFSolver::adjustPrecision(double error) {
    if (this->orbPrec[0] > 0.0) this->orbPrec[0] *= 0.5;
    this->orbPrec[0] = std::min(10.0 * error * error, this->orbPrec[0]);
    this->orbPrec[0] = std::max(this->orbPrec[0], this->orbPrec[2]);

    mrcpp::print::separator(2, '=');
    mrcpp::print::value(1, "Current precision", this->orbPrec[0], "(rel)", 5);
    mrcpp::print::separator(1, '-');
    mrcpp::print::value(2, "Orbital threshold", this->orbThrs, "(abs)", 5);
    mrcpp::print::value(2, "Property threshold", this->propThrs, "(abs)", 5);
    mrcpp::print::separator(2, '=', 2);
    return this->orbPrec[0];
}

/** @brief Get precision parameter for Helmholtz
 *
 * If the precision is NOT set explicitly (prec < 0) the current dynamic
 * precision (orbPrec[0]) will be used. If this is not set (prec < 0)
 * the current max precision (orbPrec[2]) will be used.
 */
double SCFSolver::getHelmholtzPrec() {
    double prec = this->helmPrec;
    if (prec < 0.0) prec = this->orbPrec[0];
    if (prec < 0.0) prec = this->orbPrec[2];
    if (prec < 0.0) MSG_WARN("Negative Helmholtz precision");
    return prec;
}

/** @brief Test if current errors are within the thresholds
 *
 * @param err_o: current orbital error
 * @param err_p: current property error
 *
 * A negative threshold means that it is inactive.
 */
bool SCFSolver::checkConvergence(double err_o, double err_p) const {
    bool conv_o = false;
    bool conv_p = false;
    if (std::abs(err_o) < this->orbThrs or this->orbThrs < 0.0) conv_o = true;
    if (std::abs(err_p) < this->propThrs or this->propThrs < 0.0) conv_p = true;
    return (conv_o and conv_p);
}

/** @brief Get property update
 *
 * @param vec: convergence vector
 * @param i: position in vector
 * @param absPrec: return absolute value
 *
 * Returns the difference between the i-th and (i-1)-th entry of the convergence vector.
 */
double SCFSolver::getUpdate(const std::vector<double> &vec, int i, bool absPrec) const {
    if (i < 1 or i > vec.size()) MSG_ERROR("Invalid argument");
    double E_i = vec[i - 1];
    double E_im1 = 0.0;
    if (i > 1) { E_im1 = vec[i - 2]; }
    double E_diff = E_i - E_im1;
    if (not absPrec and std::abs(E_i) > mrcpp::MachineZero) { E_diff *= 1.0 / E_i; }
    return E_diff;
}

/** @brief Pretty printing of property update
 *
 * @param name: name of property
 * @param P: current value
 * @param dP: current update
 *
 * Adds convergence status based on the property threshold.
 */
void SCFSolver::printUpdate(int plevel, const std::string &txt, double P, double dP, double thrs) const {
    int pprec = Printer::getPrecision();
    int w0 = (Printer::getWidth() - 1);
    int w1 = 25;
    int w2 = w0 / 3;
    int w3 = 8;
    int w4 = w0 - w1 - w2 - w3;

    bool done = (std::abs(dP) < thrs) or (thrs < 0.0);

    std::stringstream o_row;
    o_row << txt << std::string(w1 - txt.size(), ' ');
    o_row << std::setw(w2) << std::setprecision(2 * pprec) << std::fixed << P;
    o_row << std::setw(w4) << std::setprecision(pprec) << std::scientific << dP;
    o_row << std::setw(w3) << done;
    println(plevel, o_row.str());
}

/** @brief Pretty printing of orbitals with energies
 *
 * @param epsilon: orbital energies
 * @param Phi: orbital vector
 * @param flag: interpret epsilon as energy or norm
 */
void SCFSolver::printOrbitals(const DoubleVector &norms, const DoubleVector &errors, OrbitalVector &Phi, int flag, bool print_head) const {
    int pprec = Printer::getPrecision();
    int w0 = (Printer::getWidth() - 1);
    int w1 = 5;
    int w2 = 7;
    int w3 = 8;
    int w4 = w0 / 3;
    int w5 = 8;
    int w6 = w0 - w1 - w2 - w3 - w4 - w5;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w2) << "Spin";
    o_head << std::setw(w3) << "Nodes";
    if (flag == 0) o_head << std::setw(w4) << "F(i,i)";
    if (flag == 1) o_head << std::setw(w4) << "Norm";
    o_head << std::setw(w6) << "Residual";
    o_head << std::setw(w5) << "Done";

    if (print_head) mrcpp::print::separator(2, '=');
    if (print_head) println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    bool conv_tot = true;
    for (int i = 0; i < Phi.size(); i++) {
        bool conv_i = (errors(i) < this->orbThrs) or (this->orbThrs < 0.0);
        std::stringstream o_row;
        o_row << std::setw(w1) << i;
        o_row << std::setw(w2) << Phi[i].printSpin();
        o_row << std::setw(w3) << Phi[i].getNNodes(NUMBER::Total);
        o_row << std::setw(w4) << std::setprecision(2 * pprec) << std::fixed << norms(i);
        o_row << std::setw(w6) << std::setprecision(pprec) << std::scientific << errors(i);
        o_row << std::setw(w5) << conv_i;
        println(2, o_row.str());
        if (Phi[i].hasReal()) println(5, Phi[i].real());
        if (not conv_i) conv_tot = false;
    }
}

void SCFSolver::printResidual(double residual, bool converged) const {
    int pprec = Printer::getPrecision();
    int w0 = (Printer::getWidth() - 1);
    int w1 = 5;
    int w2 = 7;
    int w3 = 8;
    int w4 = w0 / 3;
    int w5 = 8;
    int w6 = w0 - w1 - w2 - w3 - w4 - w5;

    std::string txt = " Total residual";
    printout(1, txt << std::string(w1 + w2 + w3 - txt.size(), ' '));
    printout(1, std::setw(w4 + w6) << std::setprecision(pprec) << std::scientific << residual);
    printout(1, std::setw(w5) << converged << std::endl);
}

void SCFSolver::printConvergenceHeader(const std::string &txt) const {
    int w0 = Printer::getWidth() - 1;
    int w1 = 5;
    int w2 = 3 * w0 / 10;
    int w3 = w0 - w1 - 2 * w2;

    std::stringstream o_head;
    o_head << std::setw(w1) << "Iter";
    o_head << std::setw(w2) << "MO residual";
    o_head << std::setw(w3) << txt;
    o_head << std::setw(w2) << "Update";

    mrcpp::print::separator(0, '=');
    println(0, o_head.str());
    mrcpp::print::separator(0, '-');
}

void SCFSolver::printConvergenceRow(int i) const {
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 3 * w0 / 10;
    auto w3 = w0 - w1 - 2 * w2;

    auto residual = this->error[i];
    auto property = this->property[i];
    auto update = getUpdate(this->property, i + 1, true);

    std::stringstream o_txt;
    o_txt << std::setw(w1) << i;
    o_txt << std::setw(w2) << std::setprecision(pprec) << std::scientific << residual;
    o_txt << std::setw(w3) << std::setprecision(2 * pprec) << std::fixed << property;
    o_txt << std::setw(w2) << std::setprecision(pprec) << std::scientific << update;
    println(0, o_txt.str());
}

/** @brief Pretty printing of convergence pattern
 *
 * @param converged: convergence status
 *
 * Prints convergence in both orbitals and property.
 */
void SCFSolver::printConvergence(bool converged, const std::string &txt) const {
    auto plevel = Printer::getPrintLevel();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = (w0 - 30) / 2;
    auto w2 = (w0 - 25) / 2;

    auto nIter = this->error.size();
    if (plevel > 0) {
        printConvergenceHeader(txt);
        for (int i = 0; i < nIter; i++) printConvergenceRow(i);
    }
    mrcpp::print::separator(0, '-');
    if (converged) {
        println(0, std::string(w1, ' ') << "SCF converged in " << nIter - 1 << " iterations!");
    } else {
        println(0, std::string(w2, ' ') << "SCF did NOT converge!!!");
    }
    mrcpp::print::separator(0, '=', 2);
}

void SCFSolver::printMemory() const {
    DoubleVector mem_vec = DoubleVector::Zero(mpi::orb_size);
    mem_vec(mpi::orb_rank) = static_cast<double>(mrcpp::details::get_memory_usage());
    mpi::allreduce_vector(mem_vec, mpi::comm_orb);

    std::string mem_unit = "(kB)";
    if (mem_vec.maxCoeff() > 512.0) {
        mem_vec.array() /= 1024.0;
        mem_unit = "(MB)";
    }
    if (mem_vec.maxCoeff() > 512.0) {
        mem_vec.array() /= 1024.0;
        mem_unit = "(GB)";
    }

    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::separator(1, '-');
    mrcpp::print::header(2, "Memory usage");
    mrcpp::print::value(1, "Total memory current process", mem_vec(mpi::orb_rank), mem_unit, 2, false);
    mrcpp::print::value(2, "Maximum memory process", mem_vec.maxCoeff(), mem_unit, 2, false);
    mrcpp::print::value(2, "Minimum memory process", mem_vec.minCoeff(), mem_unit, 2, false);
    mrcpp::print::value(2, "Average memory process", mem_vec.mean(), mem_unit, 2, false);
    if (mpi::bank_size > 0 and mpi::grand_master()) {
        if (mem_unit == "(GB)") {
            mrcpp::print::value(1, "Maximum data in bank", (double)dataBank.get_maxtotalsize() / 1024, mem_unit, 2, false);
        } else {
            mrcpp::print::value(1, "Maximum data in bank", (double)dataBank.get_maxtotalsize(), "(MB)", 2, false);
        }
    }
    mrcpp::print::separator(2, '=', 2);
}

} // namespace mrchem
