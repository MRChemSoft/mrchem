#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotentialD1.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

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
    if (mpi::bank_size < 1) return;
    Timer timer;
    mpi::barrier(mpi::comm_orb);
    OrbitalVector &Phi = *this->orbitals;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) PhiBank.put_orb(i, Phi[i]);
    }
    mpi::barrier(mpi::comm_orb);
    mrcpp::print::time(4, "Setting up exchange bank", timer);
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
        if (not mpi::my_orb(phi_p)) {
            MSG_WARN("Not computing exchange contributions that are not mine");
            return out_p;
        }
        println(4, "On-the-fly exchange");
        return calcExchange(phi_p);
    } else {
        println(4, "Precomputed exchange");
        qmfunction::deep_copy(out_p, this->exchange[i]); // TODO: check if reference kan be used
        return out_p;
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

/** @brief precomputes the exchange potential
 *
 *  @param[in] phi_p input orbital
 *
 * The exchange potential is (pre)computed among the orbitals that define the operator
 */
void ExchangePotentialD1::setupInternal(double prec) {
    Timer timerT, timerS(false), t_calc(false);
    setApplyPrec(prec);
    if (this->exchange.size() != 0) MSG_ERROR("Exchange not properly cleared");

    int id_shift = 100000; // temporary shift for not colliding with existing ids
    OrbitalVector &Ex = this->exchange;
    OrbitalVector &Phi = *this->orbitals;
    BankAccount ExBank;

    // use fixed exchange_prec if set explicitly, otherwise use setup prec
    double precf = (this->exchange_prec > 0.0) ? this->exchange_prec : prec;
    // adjust precision since we sum over orbitals
    precf /= std::min(10.0, std::sqrt(1.0 * Phi.size()));

    // Initialize this->exchange and compute own diagonal elements
    Timer timerD;
    for (auto &phi_i : Phi) {
        Orbital ex_iii = phi_i.paramCopy();
        t_calc.resume();
        if (mpi::my_orb(phi_i)) calcExchange_kij(precf, phi_i, phi_i, phi_i, ex_iii);
        t_calc.stop();
        Ex.push_back(ex_iii);
    }
    mrcpp::print::time(4, "Exchange time diagonal", timerD);

    // total size (kB) of all tree saved in bank by this process so far (until read)
    int totsize = 0;
    // size allowed to be stored before going to the bank for withdrawals (5000 MB).
    int maxSize = 5000 * 1024 * mpi::bank_size / mpi::orb_size;

    int foundcount = 0;
    int N = Phi.size();
    bool use_sym = true;
    for (int j = 0; j < N; j++) {
        Orbital &phi_j = Phi[j];
        if (not mpi::my_orb(phi_j)) continue; // compute only own j
        for (int i = 0; i < N; i++) {
            if (i == j) continue;
            // compute only half of matrix (in band where i-j < size/2 %size)
            if (i > j + (N - 1) / 2 and i > j and use_sym) continue;
            if (i > j - (N + 1) / 2 and i < j and use_sym) continue;
            Orbital phi_i;
            PhiBank.get_orb(i, phi_i, 1); // fetch also own orbitals (simpler for clean up, and they are few)

            Orbital ex_jji = phi_i.paramCopy();
            Orbital ex_iij = phi_j.paramCopy();
            t_calc.resume();
            if (use_sym) {
                // compute K_iij and K_jji in one operation
                calcExchange_kij(precf, phi_i, phi_i, phi_j, ex_iij, &ex_jji);
            } else {
                // compute only own contributions (K_jji is not computed)
                calcExchange_kij(precf, phi_i, phi_i, phi_j, ex_iij);
            }
            t_calc.stop();
            double j_fac = getSpinFactor(phi_i, phi_j);
            Ex[j].add(j_fac, ex_iij);
            ex_iij.release();

            if (ex_jji.hasReal() or ex_jji.hasImag()) {
                if (not mpi::my_orb(phi_i)) {
                    // must store contribution to exchange_i in bank
                    timerS.resume();
                    totsize += ex_jji.getSizeNodes(NUMBER::Total);
                    if (ex_jji.norm() > prec) ExBank.put_orb(i + (j + 1) * N + id_shift, ex_jji);
                    timerS.stop();
                } else {
                    // must add contribution to exchange_i
                    if (ex_jji.norm() > prec) {
                        double i_fac = getSpinFactor(phi_j, phi_i);
                        Ex[i].add(i_fac, ex_jji);
                    }
                }
            }
            ex_jji.release();

            // Periodically the bank is visited and all contributions to own exchange that are stored there
            // and ready to use are fetched, added to exchange (and deleted from bank).
            // This must not be done too often for not waisting time, but often enough so that the bank does not fill
            // too much.
            if (totsize > maxSize) { // rough estimate of how full the bank is
                timerS.resume();
                totsize = 0; // reset counter
                Orbital ex_rcv;
                for (int jj = 0; jj < N; jj++) {
                    if (not mpi::my_orb(Phi[jj])) continue; // compute only own j
                    for (int ii = 0; ii < N; ii++) {
                        if (ii == jj) continue;
                        if (mpi::my_orb(Phi[ii])) continue; // fetch only other's i
                        // Fetch other half of matrix (in band where i-j > size/2 %size)
                        if ((ii > jj + (N - 1) / 2 and ii > jj) or (ii > jj - (N + 1) / 2 and ii < jj)) {
                            double j_fac = getSpinFactor(Phi[ii], Phi[jj]);
                            int found = ExBank.get_orb_del(jj + (ii + 1) * N + id_shift, ex_rcv);
                            foundcount += found;
                            if (found) Ex[jj].add(j_fac, ex_rcv);
                        }
                    }
                    Ex[jj].crop(prec);
                }
                printout(4, " fetched " << foundcount << " Exchange contributions from bank");
                printout(4, (int)((float)timerS.elapsed() * 1000) << " ms\n");
                timerS.stop();
            }
        }
    }

    mrcpp::print::time(3, "Time exchanges compute", timerT);
    mpi::barrier(mpi::comm_orb);
    mrcpp::print::time(3, "Time exchanges all mpi finished", timerT);

    timerS.resume();
    for (int j = 0; j < N and use_sym; j++) {  // If symmetri is not used, there is nothing to fetch in bank
        if (not mpi::my_orb(Phi[j])) continue; // fetch only own j
        Orbital ex_rcv;
        for (int i = 0; i < N; i++) {
            if (i == j) continue;
            if (mpi::my_orb(Phi[i])) continue; // fetch only other's i
            // Fetch other half of matrix (in band where i-j > size/2 %size)
            if ((i > j + (N - 1) / 2 and i > j) or (i > j - (N + 1) / 2 and i < j)) {
                double j_fac = getSpinFactor(Phi[i], Phi[j]);
                int found = ExBank.get_orb_del(j + (i + 1) * N + id_shift, ex_rcv);
                foundcount += found;
                if (found) Ex[j].add(j_fac, ex_rcv);
            }
        }
        Ex[j].crop(prec);
    }
    timerS.stop();
    println(3, " fetched in total " << foundcount << " Exchange contributions from bank");
    mrcpp::print::time(3, "Time send/rcv exchanges", timerS);
    mrcpp::print::time(3, "Time calculate exchanges", t_calc);

    auto n = orbital::get_n_nodes(Ex);
    auto m = orbital::get_size_nodes(Ex);
    auto t = timerT.elapsed();
    mrcpp::print::tree(2, "Hartree-Fock exchange", n, m, t);
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

    QMFunctionVector func_vec;
    std::vector<ComplexDouble> coef_vec;
    for (int i = 0; i < Phi.size(); i++) {
        Orbital &phi_i = Phi[i];
        if (not mpi::my_orb(phi_i)) PhiBank.get_orb(i, phi_i, 1);

        double spin_fac = getSpinFactor(phi_i, phi_p);
        if (std::abs(spin_fac) >= mrcpp::MachineZero) {
            Orbital ex_iip = phi_p.paramCopy();
            calcExchange_kij(precf, phi_i, phi_i, phi_p, ex_iip);
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
            func_vec.push_back(ex_iip);
        }

        if (not mpi::my_orb(phi_i)) phi_i.free(NUMBER::Total);
    }

    // compute ex_p = sum_i c_i*ex_iip
    Orbital ex_p = phi_p.paramCopy();
    Eigen::Map<ComplexVector> coefs(coef_vec.data(), coef_vec.size());
    qmfunction::linear_combination(ex_p, coefs, func_vec, prec);
    print_utils::qmfunction(3, "Applied exchange", ex_p, timer);
    return ex_p;
}

} // namespace mrchem
