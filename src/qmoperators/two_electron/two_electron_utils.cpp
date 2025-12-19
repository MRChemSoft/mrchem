#include "two_electron_utils.h"
#include "qmfunctions/density_utils.h"
#include <MRCPP/Parallel>
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

using mrcpp::Printer;
using mrcpp::Timer;
extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

namespace mrchem {

/** @brief precomputes all the <ik|1/r12|jl> 2-electron integrals
 *
 *  @param[in] Phi input orbitals
 *
 * The ij/r12 potentials are first (pre)computed for all the orbitals j<i, then multiplied by kl
 */
Eigen::Tensor<std::complex<double>, 4> calc_2elintegrals(double prec, OrbitalVector &Phi) {
    mrcpp::print::header(3, "Computing two electron integrals");
    Timer t_tot, t_snd(false), t_rho(false), t_ovr(false), t_calc(false), t_add(false), t_get(false);

    mrcpp::BankAccount PhiBank; // to put the orbitals
    mrcpp::BankAccount VijBank; // to put the ij/r12 potentials

    int N = Phi.size();
    double precf = prec / std::sqrt(1 * Phi.size());
    double prec_m2 = precf / 100; // second multiplication

    for (int i = 0; i < Phi.size(); i++) {
        if (mrcpp::mpi::my_func(i)) PhiBank.put_func(i, Phi[i]);
    }

    std::vector<mrcpp::CompFunction<3>> Vij_vec;
    if (mrcpp::mpi::bank_size <= 0) {
        Vij_vec.resize(N*(N+1)/2);
    }

    // first we compute density to be used so we know were the orbitals are. This
    // is to avoid computing the potential at positions were there are no orbitals.
    t_rho.resume();
    Density rho;
    density::compute(prec, rho, Phi, DensityType::Total);
    t_rho.stop();

    // Initialize and compute own diagonal elements
    Timer t_diag;
    int i = -1;
    for (auto &phi_i : Phi) {
        i++;
        if (!mrcpp::mpi::my_func(i)) continue;
        Orbital Vij = phi_i.paramCopy(true);
        t_calc.resume();
        ij_r12(precf, rho, phi_i,  phi_i, Vij);
        t_calc.stop();
        if (mrcpp::mpi::bank_size > 0) {
            // store Vij
            t_snd.resume();
            VijBank.put_func(i*(i+1)/2+i, Vij);
            // if (Vij.norm() > prec)VijBank.put_func(i*(i+1)/2+i, Vij);
            t_snd.stop();
        } else {
            Vij_vec[i*(i+1)/2+i] = Vij;
        }
        Vij.free();
    }
    t_diag.stop();

    Timer t_offd;
    // We divide all the contributions into a fixed number of tasks.
    // all "j" orbitals are fetched and stored, and used together with one "i" orbital

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

    // store Int(phi_i^dag*phi_j/|r-r'|) potentials in Bank
    mrcpp::TaskManager tasksMaster(ntasks);
    while (true) {
        task = tasksMaster.next_task();
        if (task < 0) break;
        // we fetch all required i (but only one j at a time)
        std::vector<Orbital> iorb_vec;
        int i0 = -1;
        for (int i = 0; i < itasks[task].size(); i++) {
            int iorb = itasks[task][i];
            i0 = iorb;
            Orbital phi_i;
            t_get.resume();
            if (mrcpp::mpi::bank_size > 0) {
                PhiBank.get_func(iorb, phi_i, 1); // fetch also own orbitals (simpler for clean up, and they are few)
                iorb_vec.push_back(phi_i);
            } else {
                iorb_vec.push_back(Phi[iorb]);
            }
            t_get.stop();
        }

        for (int j = 0; j < jtasks[task].size(); j++) {
            int jorb = jtasks[task][j];
            Orbital phi_j;
            t_get.resume();
            if (mrcpp::mpi::bank_size > 0) {
                PhiBank.get_func(jorb, phi_j, 1);
            } else
                phi_j = Phi[jorb];
            t_get.stop();

            for (int i = 0; i < iorb_vec.size(); i++) {
                int iorb = itasks[task][i];
                Orbital &phi_i = iorb_vec[i];
                Orbital Vij = phi_i.paramCopy(true);
                t_calc.resume();
                ij_r12(precf, rho, phi_i, phi_j, Vij);
                t_calc.stop();
                int ij = std::max(iorb,jorb)*(std::max(iorb,jorb)+1)/2 + std::min(iorb,jorb);
                if (mrcpp::mpi::bank_size > 0) {
                    // store Vij
                    t_snd.resume();
                    VijBank.put_func(ij, Vij);
                    // if (Vij.norm() > prec)VijBank.put_func(ij, Vij);
                    t_snd.stop();
                } else {
                    Vij_vec[ij] = Vij;
                }
                Vij.free();
            }
        }
    }
    t_offd.stop();
    Eigen::Tensor<std::complex<double>, 4> klVij(N,N,N,N);
    // For each Vij, compute the overlap matrix <k| lVij>
    long long totSizeVij = 0;
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j <= i; j ++) {
            Orbital Vij;
            t_get.resume();
            VijBank.get_func(i*(i+1)/2+j, Vij, 1); // note that it will wait if not yet ready
            totSizeVij += Vij.getSizeNodes(); // in kB
            t_get.stop();

            //compute l*Vij for own orbitals only
            OrbitalVector l_Vij;
            for (auto &phi_l : Phi) {
                Orbital lVij = phi_l.paramCopy(true);
                if (mrcpp::mpi::my_func(phi_l)) {
                    mrcpp::multiply(lVij, phi_l, Vij, prec_m2, true, true);
                }
                l_Vij.push_back(lVij);
            }
            Vij.free();
            t_ovr.resume();
            ComplexMatrix kl_Vij = calc_overlap_matrix(Phi,l_Vij);
            t_ovr.stop();

            for (int k = 0; k < N; k ++) {
                for (int l = 0; l < N; l ++) {
                    klVij(i,j,k,l) = kl_Vij(k,l);
                    klVij(j,i,k,l) = kl_Vij(k,l);
                }
            }
        }
    }
    t_tot.stop();

    mrcpp::print::time(3, "Time making density", t_rho);
    mrcpp::print::time(3, "Time receiving orbitals", t_get);
    mrcpp::print::time(3, "Time sending potentials", t_snd);
    mrcpp::print::time(3, "Time computing integrals", t_calc);
    mrcpp::print::separator(3, '-');
    mrcpp::print::time(3, "Time diagonal terms", t_diag);
    mrcpp::print::time(3, "Time off-diagonal terms", t_offd);
    mrcpp::print::time(3, "Time making overlap matrix", t_ovr);
    mrcpp::print::time(3, "Time making two el integrals", t_tot);
    mrcpp::print::value(3,"Total size all ij/r12 potentials", totSizeVij/1024.0, "(MB)",0,false);
    mrcpp::print::separator(3, '-');

    auto t = t_tot.elapsed() / N;
    return klVij;
}

/** @brief computes Int(phi_i^dag*phi_j/|r-r'|)
 *
 *  \param[in] phi_i orbital to be conjugated and multiplied by phi_j
 *  \param[in] phi_j orbital to be multiplied by phi_i^dag
 *  \param[in] rho is normally the density. It is used to screen the final output
 *  \param[out] V_ij is the resulting potential function
  *
 * Computes the product of complex conjugate of phi_i and phi_j,
 * then applies the Poisson operator. The result is given in V_ij.
 */
void ij_r12(double prec, Orbital rho, Orbital phi_i, Orbital phi_j, Orbital &V_ij) {
    Timer timer_tot;
    mrcpp::PoissonOperator P(*MRA, prec);

    // set precisions
    double prec_m1 = prec / 10;  // first multiplication
    double prec_p = prec * 10;   // Poisson application

    // compute rho_ij = phi_i^dagger * phi_j
    // if the product is smaller than the target precision,
    // the result is expected to be negligible
    Timer timer_ij;
    Orbital rho_ij = phi_i.paramCopy(true);
    mrcpp::multiply(rho_ij, phi_i, phi_j, prec_m1, true, true, true);
    timer_ij.stop();
    if (rho_ij.norm() < prec) return;

    auto N_i = phi_i.getNNodes();
    auto N_j = phi_j.getNNodes();
    auto N_ij = rho_ij.getNNodes();
    auto norm_ij = rho_ij.norm();
    // For now we assume all phi are complex or all ar real.

    bool RealOrbitals = phi_i.isreal();

    // prepare vector used to steer precision of Poisson application
    mrcpp::FunctionTreeVector<3, double> phi_opt_vec_real;
    mrcpp::FunctionTreeVector<3, ComplexDouble> phi_opt_vec_cplx;
    if (RealOrbitals) {
        if (rho.isreal() and rho.getNNodes() > 0) phi_opt_vec_real.emplace_back(1.0, rho.CompD[0]);
        if (phi_j.isreal() and  phi_j.getNNodes() > 0) phi_opt_vec_real.emplace_back(1.0, phi_j.CompD[0]);
        if (phi_i.isreal() and  &phi_i != &phi_j and phi_i.getNNodes() > 0) phi_opt_vec_real.emplace_back(1.0, phi_i.CompD[0]);
    } else {
        if (rho.iscomplex()) phi_opt_vec_cplx.emplace_back(1.0, rho.CompC[0]);
        if (phi_j.iscomplex()) phi_opt_vec_cplx.emplace_back(1.0, phi_j.CompC[0]);
        if (phi_i.iscomplex() and &phi_i != &phi_j) phi_opt_vec_cplx.emplace_back(1.0, phi_i.CompC[0]);
    }
    // compute V_ij = P[rho_ij]
    Timer timer_p;
    if (RealOrbitals) {
        mrcpp::apply(prec_p, *V_ij.CompD[0], P, *rho_ij.CompD[0], phi_opt_vec_real, -1, true);
    } else {
        mrcpp::apply(prec_p, *V_ij.CompC[0], P, *rho_ij.CompC[0], phi_opt_vec_cplx, -1, true);
    }
    rho_ij.free();
    timer_p.stop();
    auto N_p = V_ij.getNNodes();
    auto norm_p = V_ij.norm();

    println(5,
            " time " << (int)((float)timer_tot.elapsed() * 1000) << " ms "
            << " mult1:" << (int)((float)timer_ij.elapsed() * 1000) << " Pot:" << (int)((float)timer_p.elapsed() * 1000) << " Nnodes: " << N_i << " " << N_j << " " << N_ij << " " << N_p << " " << " norms " << norm_ij<<" "
                     << norm_p );
}
}
