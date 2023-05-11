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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotentialD1.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;
using QMOperator_p = std::shared_ptr<mrchem::QMOperator>;

namespace mrchem {

/** @brief constructor
 *
 * @param[in] P Poisson operator (does not take ownership)
 * @param[in] Phi vector of orbitals which define the exchange operator
 * @param[in] prec screening precision for exchange construction
 */
ExchangePotentialD1::ExchangePotentialD1(PoissonOperator_p P, OrbitalVector_p Phi, double prec)
        : ExchangePotential(P, Phi, prec) {}

/** @brief Save all orbitals in Bank, so that they can be accessed asynchronously */
void ExchangePotentialD1::setupBank() {
    if (mrcpp::mpi::bank_size < 1) return;
    Timer timer;
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    OrbitalVector &Phi = *this->orbitals;
    for (int i = 0; i < Phi.size(); i++) {
        if (mrcpp::mpi::my_orb(i)) PhiBank.put_func(i, Phi[i]);
    }
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    mrcpp::print::time(3, "Setting up exchange bank", timer);
}

/** @brief Clears rbital bank accounts.
 *
 */
void ExchangePotentialD1::clearBank() {
    PhiBank.clear();
}

/** @brief Test if a given contribution has been precomputed
 *
 * @param[in] phi_p orbital for which the check is performed
 *
 * If the given contribution has been precomputed, it is simply copied,
 * without additional recalculation.
 */
int ExchangePotentialD1::testInternal(Orbital phi_p) const {
    const OrbitalVector &Phi = *this->orbitals;
    const OrbitalVector &Kphi = this->exchange;

    int out = -1;
    if (Kphi.size() == Phi.size()) {
        for (int i = 0; i < Phi.size(); i++) {
            if (&Phi[i].real() == &phi_p.real() and &Phi[i].imag() == &phi_p.imag()) {
                out = i;
                break;
            }
        }
    }
    return out;
}

/** @brief Applies operator potential
 *
 *  @param[in] inp input orbital
 *
 * The exchange potential is applied to the given orbital. Checks first if this
 * particular exchange contribution has been precomputed.
 */
Orbital ExchangePotentialD1::apply(Orbital phi_p) {
    Orbital out_p = phi_p.paramCopy();
    if (this->apply_prec < 0.0) {
        MSG_ERROR("Uninitialized operator");
        return out_p;
    }
    int i = testInternal(phi_p);
    if (i < 0) {
        if (not mrcpp::mpi::my_orb(phi_p)) {
            MSG_WARN("Not computing exchange contributions that are not mine");
            return out_p;
        }
        return calcExchange(phi_p);
    } else {
        return this->exchange[i];
    }
}

/** @brief Applies the adjoint of the operator
 *  \param[in] inp input orbital
 *
 *  Self-adjoint operator
 */
Orbital ExchangePotentialD1::dagger(Orbital phi_p) {
    return apply(phi_p);
}

QMOperatorVector ExchangePotentialD1::apply(QMOperator_p &O) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief precomputes the exchange potential
 *
 *  @param[in] phi_p input orbital
 *
 * The exchange potential is (pre)computed among the orbitals that define the operator
 */
void ExchangePotentialD1::setupInternal(double prec) {
    Timer t_tot, t_snd(false), t_orb(false), t_calc(false), t_add(false), t_get(false), t_wait(false);
    setApplyPrec(prec);
    if (this->exchange.size() != 0) MSG_ERROR("Exchange not properly cleared");

    OrbitalVector &Ex = this->exchange;
    OrbitalVector &Phi = *this->orbitals;
    mrcpp::BankAccount ExBank;
    int N = Phi.size();
    // use fixed exchange_prec if set explicitly, otherwise use setup prec
    double precf = (this->exchange_prec > 0.0) ? this->exchange_prec : prec;
    prec = mrcpp::mpi::numerically_exact ? -1.0 : prec;
    precf /= std::sqrt(1 * Phi.size());
    // Initialize this->exchange and compute own diagonal elements
    Timer t_diag;
    int i = 0;
    for (auto &phi_i : Phi) {
        Orbital ex_iii(phi_i.spin(), phi_i.occ(), phi_i.getRank());
        t_calc.resume();
        if (mrcpp::mpi::my_orb(i)) calcExchange_kij(precf, phi_i, phi_i, phi_i, ex_iii);
        t_calc.stop();
        Ex.push_back(ex_iii);
        i++;
    }
    t_diag.stop();

    Timer t_offd;
    // We divide all the exchange contributions into a fixed number of tasks.
    // all "j" orbitals are fetched and stored, and used together with one "i" orbital
    // At the end of a task, the Exchange for j is summed up with the contributions
    // stored in the bank, and then stored.
    // When all tasks are finished, each MPI fetch and sum the contributions to its
    // own exchange.

    // make a set of tasks
    // We use symmetry: each pair (i,j) must be used once only. Only j<i
    // Divide into square blocks, with the diagonal blocks taken at the end (because they are faster to compute)
    int block_size; // NB: block_size*block_size intermediate exchange results are stored temporarily
    block_size = std::min(16, std::max(2, static_cast<int>(std::sqrt(N * N / (14 * mrcpp::mpi::wrk_size)))));

    int iblocks = (N + block_size - 1) / block_size;
    int ntasksmax = ((iblocks - 1) * iblocks) / 2 + iblocks * (block_size * (block_size - 1) / 2);
    std::vector<std::vector<int>> itasks(ntasksmax); // the i values (orbitals) of each block
    std::vector<std::vector<int>> jtasks(ntasksmax); // the j values (orbitals) of each block

    int task = 0;
    // make a path for tasks that follow diagonals, in order to maximize the spread of orbitals treated
    for (int ib = 0; ib < (iblocks + 1) / 2; ib++) {
        if (task >= (iblocks * (iblocks - 1) / 2)) break;
        for (int ij = 0; ij < iblocks; ij++) {
            int j0 = ij;
            int i0 = ib + ij + 1;
            if (i0 * block_size >= N) {
                // continue in "symmetrical" part
                j0 = ij + ib + 1 - iblocks;
                i0 = ij;
            }
            for (int jj = 0; jj < block_size; jj++) {
                int jjj = j0 * block_size + jj;
                if ((i0 + j0) % 2 != 0) jjj = j0 * block_size + (block_size - 1 - jj); // reversed order
                if (jjj >= N) continue;
                if ((i0 + j0) % 2 == 0)
                    jtasks[task].push_back(jjj);
                else
                    itasks[task].push_back(jjj);
            }
            for (int ii = 0; ii < block_size; ii++) {
                int iii = i0 * block_size + ii;
                if ((i0 + j0) % 2 == 0) iii = i0 * block_size + (block_size - 1 - ii); // reversed order
                if (iii >= N) continue;

                if ((i0 + j0) % 2 == 0)
                    itasks[task].push_back(iii);
                else
                    jtasks[task].push_back(iii);
            }
            task++;
            if (task >= (iblocks * (iblocks - 1) / 2)) break;
        }
    }

    // add diagonal blocks:
    // we make those tasks smaller (1x1 blocks), in order to minimize the time waiting for the last task.
    // NB: only include j<i within these blocks
    for (int i = 0; i < N; i += block_size) {
        int j = i;
        for (int jj = j; jj < j + block_size and jj < N; jj++) {
            for (int ii = i; ii < i + block_size and ii < N; ii++) {
                if (ii <= jj) continue; // only jj<ii is computed
                itasks[task].push_back(ii);
                jtasks[task].push_back(jj);
                task++;
            }
        }
    }
    assert(task <= ntasksmax);
    int ntasks = task;

    mrcpp::TaskManager tasksMaster(ntasks);
    while (true) {
        task = tasksMaster.next_task();
        if (task < 0) break;
        // we fetch all required i (but only one j at a time)
        OrbitalVector iorb_vec;
        int i0 = -1;
        for (int i = 0; i < itasks[task].size(); i++) {
            int iorb = itasks[task][i];
            i0 = iorb;
            Orbital phi_i;

            t_orb.resume();
            if (mrcpp::mpi::bank_size > 0) {
                PhiBank.get_func(iorb, phi_i, 1); // fetch also own orbitals (simpler for clean up, and they are few)
                iorb_vec.push_back(phi_i);
            } else {
                iorb_vec.push_back(Phi[iorb]);
            }
            t_orb.stop();
        }

        for (int j = 0; j < jtasks[task].size(); j++) {
            int jorb = jtasks[task][j];
            Orbital phi_j;
            t_orb.resume();
            if (mrcpp::mpi::bank_size > 0) {
                PhiBank.get_func(jorb, phi_j, 1);
            } else
                phi_j = Phi[jorb];
            t_orb.stop();
            std::vector<mrcpp::ComplexFunction> iijfunc_vec;
            ComplexVector coef_vec(N);
            for (int i = 0; i < iorb_vec.size(); i++) {
                int iorb = itasks[task][i];
                Orbital &phi_i = iorb_vec[i];
                Orbital ex_jji = phi_i.paramCopy();
                Orbital ex_iij = phi_j.paramCopy();

                // compute K_iij and K_jji in one operation
                double j_fac = getSpinFactor(phi_i, phi_j);
                if (std::abs(j_fac) < mrcpp::MachineZero) continue;
                t_calc.resume();
                calcExchange_kij(precf, phi_i, phi_i, phi_j, ex_iij, &ex_jji);
                t_calc.stop();
                if (ex_iij.norm() > prec) coef_vec[iijfunc_vec.size()] = j_fac;
                t_snd.resume();
                if (mrcpp::mpi::bank_size > 0) {
                    // store ex_jji
                    if (ex_iij.norm() > prec) iijfunc_vec.push_back(ex_iij);
                    if (ex_jji.norm() > prec) ExBank.put_func(iorb + jorb * N, ex_jji);
                    if (ex_jji.norm() > prec) tasksMaster.put_readytask(iorb, jorb);
                } else {
                    Ex[iorb].add(j_fac, ex_jji);
                    Ex[jorb].add(j_fac, ex_iij);
                }
                ex_jji.free(NUMBER::Total);
                t_snd.stop();
            }
            Timer timerx;
            // fetch ready contributions to ex_j from others
            std::vector<int> iVec = tasksMaster.get_readytask(jorb, 1);
            int lastsize = iVec.size();
            for (int iorb : iVec) {
                t_get.resume();
                Orbital ex_rcv;
                int found = ExBank.get_func_del(jorb + iorb * N, ex_rcv);
                t_get.stop();
                if (not found) MSG_ERROR("Exchange not found");
                double j_fac = getSpinFactor(ex_rcv, phi_j);
                coef_vec[iijfunc_vec.size()] = j_fac;
                iijfunc_vec.push_back(ex_rcv);
            }
            // add all contributions to ex_j,
            if (mrcpp::mpi::bank_size > 0 and iijfunc_vec.size() > 0) {
                Orbital ex_j = phi_j.paramCopy();
                t_add.resume();
                mrcpp::cplxfunc::linear_combination(ex_j, coef_vec, iijfunc_vec, prec);
                t_add.stop();
                // ex_j is sent to Bank
                if (ex_j.hasReal() or ex_j.hasImag()) {
                    auto tT = timerx.elapsed();
                    t_snd.resume();
                    ex_j.crop(prec);
                    if (ex_j.norm() > prec) {
                        if (i0 < 0) MSG_ERROR("no exchange contributions?");
                        ExBank.put_func(jorb + (i0 + N) * N, ex_j);
                        tasksMaster.put_readytask(jorb, i0 + N);
                    }
                    t_snd.stop();
                    ex_j.free(NUMBER::Total);
                } else if (iijfunc_vec.size() > 0) {
                    MSG_ERROR("Exchange exists but has no real and no Imag parts");
                }
                for (int jj = 0; jj < iijfunc_vec.size(); jj++) iijfunc_vec[jj].free(NUMBER::Total);
            }
        }
    }

    // wait until all exchanges pieces are computed and stored in Bank
    t_wait.resume();
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    t_wait.stop();

    for (int j = 0; j < N; j++) {
        if (not mrcpp::mpi::my_orb(j) or mrcpp::mpi::bank_size == 0) continue; // fetch only own j
        std::vector<int> iVec = tasksMaster.get_readytask(j, 1);
        std::vector<mrcpp::ComplexFunction> iijfunc_vec;
        ComplexVector coef_vec(N);
        int tot = 0;
        int totmax = 2 * block_size;
        for (int i : iVec) {
            if (i < 0) continue;
            t_get.resume();
            Orbital ex_rcv;
            int found = ExBank.get_func_del(j + i * N, ex_rcv);
            t_get.stop();
            double j_fac = getSpinFactor(ex_rcv, Phi[j]);
            coef_vec[iijfunc_vec.size()] = j_fac;
            iijfunc_vec.push_back(ex_rcv);
            if (not found) MSG_ERROR("My Exchange not found in Bank");
            tot++;
            if (tot >= totmax) {
                // we sum the contributions so far before fetching new ones
                t_add.resume();
                auto tmp_j = Ex[j].paramCopy();
                mrcpp::cplxfunc::linear_combination(tmp_j, coef_vec, iijfunc_vec, prec);
                Ex[j].add(1.0, tmp_j);
                tmp_j.free(NUMBER::Total);
                for (int jj = 0; jj < iijfunc_vec.size(); jj++) iijfunc_vec[jj].free(NUMBER::Total);
                iijfunc_vec.clear();
                Ex[j].crop(prec);
                t_add.stop();
                tot = 0;
            }
        }
        if (iijfunc_vec.size() > 0) {
            t_add.resume();
            auto tmp_j = Ex[j].paramCopy();
            mrcpp::cplxfunc::linear_combination(tmp_j, coef_vec, iijfunc_vec, prec);
            Ex[j].add(1.0, tmp_j);
            tmp_j.free(NUMBER::Total);
            for (int jj = 0; jj < iijfunc_vec.size(); jj++) iijfunc_vec[jj].free(NUMBER::Total);
            Ex[j].crop(prec);
            t_add.stop();
        }
    }
    t_offd.stop();
    mrcpp::print::time(3, "Time receiving orbitals", t_orb);
    mrcpp::print::time(3, "Time receiving exchanges", t_get);
    mrcpp::print::time(3, "Time sending exchanges", t_snd);
    mrcpp::print::time(3, "Time adding exchanges", t_add);
    mrcpp::print::time(3, "Time waiting for others", t_wait);
    mrcpp::print::time(3, "Time computing exchanges", t_calc);
    mrcpp::print::separator(3, '-');
    mrcpp::print::time(3, "Time diagonal terms", t_diag);
    mrcpp::print::time(3, "Time off-diagonal terms", t_offd);
    mrcpp::print::separator(3, '-');

    auto t = t_tot.elapsed() / N;
    auto n = orbital::get_n_nodes(this->exchange, true);
    auto m = orbital::get_size_nodes(this->exchange, true);
    mrcpp::print::tree(3, "Average exchange term", n, m, t);
}

/** @brief Computes the exchange potential on the fly
 *
 *  \param[in] phi_p input orbital
 *
 * The exchange potential is computed and applied on the fly to the given orbital.
 * The orbitals must have been previously stored in bank.
 */
Orbital ExchangePotentialD1::calcExchange(Orbital phi_p) {
    Timer timer;
    OrbitalVector &Phi = *this->orbitals;

    double prec = this->apply_prec;
    // use fixed exchange_prec if set explicitly, otherwise use setup prec
    double precf = (this->exchange_prec > 0.0) ? this->exchange_prec : prec;
    // adjust precision since we sum over orbitals
    precf /= std::min(10.0, std::sqrt(1.0 * Phi.size()));

    std::vector<mrcpp::ComplexFunction> func_vec;
    std::vector<ComplexDouble> coef_vec;
    for (int i = 0; i < Phi.size(); i++) {
        Orbital &phi_i = Phi[i];
        if (not mrcpp::mpi::my_orb(i)) PhiBank.get_func(i, phi_i, 1);

        double spin_fac = getSpinFactor(phi_i, phi_p);
        if (std::abs(spin_fac) >= mrcpp::MachineZero) {
            Orbital ex_iip = phi_p.paramCopy();
            calcExchange_kij(precf, phi_i, phi_i, phi_p, ex_iip);
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
            func_vec.push_back(ex_iip);
        }

        if (not mrcpp::mpi::my_orb(i)) phi_i.free(NUMBER::Total);
    }

    // compute ex_p = sum_i c_i*ex_iip
    Orbital ex_p = phi_p.paramCopy();
    Eigen::Map<ComplexVector> coefs(coef_vec.data(), coef_vec.size());
    mrcpp::cplxfunc::linear_combination(ex_p, coefs, func_vec, prec);
    print_utils::qmfunction(4, "Applied exchange", ex_p, timer);
    return ex_p;
}

} // namespace mrchem
