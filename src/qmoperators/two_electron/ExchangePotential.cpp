#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotential.h"
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
ExchangePotential::ExchangePotential(PoissonOperator_p P, OrbitalVector_p Phi, double prec)
        : exchange_prec(prec)
        , orbitals(Phi)
        , poisson(P) {}

/** @brief Perform a unitary transformation among the precomputed exchange contributions
 *
 * @param[in] U unitary matrix defining the rotation
 */
void ExchangePotential::rotate(const ComplexMatrix &U) {
    if (this->exchange.size() == 0) return;
    this->exchange = orbital::rotate(this->exchange, U, this->apply_prec);

    // NOTE: The following MPI point is currently NOT implemented!
    //
    // the last parameter, 1, means MPI will send only one orbital at a time
    // (because Exchange orbitals can be large for large molecules).
    // OrbitalAdder add(this->apply_prec, this->max_scale, 1);
    // add.rotate(this->exchange, U);
}

/** @brief determines the exchange factor to be used in the calculation of the exact exchange
 *
 * @param [in] phi_i orbital defining the K operator
 * @param [in] phi_j orbital to which K is applied
 *
 * The factor is computed in terms of the occupancy of the two orbitals and in terms of the spin
 * 0.5 factors are used in order to preserve occupancy of the set of doubly occupied orbitals
 * this-> is the orbital defining the operator whereas the input orbital (orb) is the one
 * the operator is applied to
 *
 * Occupancy: Single/Double
 * Spin: alpha/beta
 *
 * K (this->) | orb (input) | factor
 * alpha      | alpha       | 1.0
 * alpha      | beta        | 0.0
 * alpha      | double      | 0.5
 * -------------------------------
 * beta       | alpha       | 0.0
 * beta       | beta        | 1.0
 * beta       | double      | 0.5
 * -------------------------------
 * double     | alpha       | 1.0
 * double     | beta        | 1.0
 * double     | double      | 1.0
 *
 */
double ExchangePotential::getSpinFactor(Orbital phi_i, Orbital phi_j) const {
    double out = 0.0;
    if (phi_i.spin() == SPIN::Paired)
        out = 1.0;
    else if (phi_j.spin() == SPIN::Paired)
        out = 0.5;
    else if (phi_i.spin() == phi_j.spin())
        out = 1.0;
    return out;
}

/** @brief Prepare operator for application
 *
 * @param[in] prec reqested precision
 *
 * This will save the defining orbitals in the orbital bank
 * and optionally compute the internal constributions.
 */
void ExchangePotential::setup(double prec) {
    if (mpi::world_size > 1 and mpi::bank_size < 1) MSG_ABORT("MPI bank required!");
    setApplyPrec(prec);
    setupBank();
    if (this->pre_compute) setupInternal(prec);
}

/** @brief Clears the Exchange Operator
 *
 *  Deletes the precomputed exchange contributions and
 *  clears the orbital bank.
 */
void ExchangePotential::clear() {
    clearInternal();
    clearBank();
    clearApplyPrec();
}

/** @brief computes phi_k*Int(phi_i^dag*phi_j/|r-r'|)
 *
 *  \param[in] phi_k orbital to be multiplied after application of Poisson operator
 *  \param[in] phi_i orbital to be conjugated and multiplied by phi_j
 *  \param[in] phi_j orbital to be multiplied by phi_i^dag
 *  \param[out] out_kij result
 *  \param[out] out_jji (optional), result where phi_k is replaced by phi_j (i.e. phi_k not used, and phi_j used
 * twice)
 *
 * Computes the product of complex conjugate of phi_i and phi_j,
 * then applies the Poisson operator, and multiplies the result
 * by phi_k (and optionally by phi_j). The result is given in phi_out.
 */
void ExchangePotential::calcExchange_kij(double prec,
                                         Orbital phi_k,
                                         Orbital phi_i,
                                         Orbital phi_j,
                                         Orbital &out_kij,
                                         Orbital *out_jji) {
    Timer timer_tot;
    mrcpp::PoissonOperator &P = *this->poisson;

    // set precisions
    double prec_m1 = prec / 10;  // first multiplication
    double prec_p = prec * 10;   // Poisson application
    double prec_m2 = prec / 100; // second multiplication

    // compute rho_ij = phi_i^dagger * phi_j
    // if the product is smaller than the target precision,
    // the result is expected to be negligible
    Timer timer_ij;
    Orbital rho_ij = phi_i.paramCopy();
    qmfunction::multiply(rho_ij, phi_i.dagger(), phi_j, prec_m1, true, true);
    timer_ij.stop();
    if (rho_ij.norm() < prec) return;

    auto N_i = phi_i.getNNodes(NUMBER::Total);
    auto N_j = phi_j.getNNodes(NUMBER::Total);
    auto N_ij = rho_ij.getNNodes(NUMBER::Total);
    auto norm_ij = rho_ij.norm();

    // prepare vector used to steer precision of Poisson application
    mrcpp::FunctionTreeVector<3> phi_opt_vec;
    if (phi_k.hasReal()) phi_opt_vec.push_back(std::make_tuple(1.0, &phi_k.real()));
    if (phi_k.hasImag()) phi_opt_vec.push_back(std::make_tuple(1.0, &phi_k.imag()));
    if (out_jji != nullptr) {
        if (phi_j.hasReal()) phi_opt_vec.push_back(std::make_tuple(1.0, &phi_j.real()));
        if (phi_j.hasImag()) phi_opt_vec.push_back(std::make_tuple(1.0, &phi_j.imag()));
    }

    // compute V_ij = P[rho_ij]
    Timer timer_p;
    Orbital V_ij = rho_ij.paramCopy();
    if (rho_ij.hasReal()) {
        V_ij.alloc(NUMBER::Real);
        mrcpp::apply(prec_p, V_ij.real(), P, rho_ij.real(), phi_opt_vec, -1, true);
    }
    if (rho_ij.hasImag()) {
        V_ij.alloc(NUMBER::Imag);
        mrcpp::apply(prec_p, V_ij.imag(), P, rho_ij.imag(), phi_opt_vec, -1, true);
    }
    rho_ij.release();
    timer_p.stop();
    auto N_p = V_ij.getNNodes(NUMBER::Total);
    auto norm_p = V_ij.norm();

    // compute out_kij = phi_k * V_ij
    Timer timer_kij;
    qmfunction::multiply(out_kij, phi_k, V_ij, prec_m2, true, true);
    auto N_kij = out_kij.getNNodes(NUMBER::Total);
    auto norm_kij = out_kij.norm();
    timer_kij.stop();

    // compute out_jji = phi_j * V_ji = phi_j * V_ij^dagger
    Timer timer_jji;
    auto N_jji = 0;
    auto norm_jji = 0.0;
    if (out_jji != nullptr) {
        qmfunction::multiply(*out_jji, phi_j, V_ij.dagger(), prec_m2, true, true);
        N_jji = out_jji->getNNodes(NUMBER::Total);
        norm_jji = out_jji->norm();
    }
    timer_jji.stop();

    println(4,
            " time " << (int)((float)timer_tot.elapsed() * 1000) << " ms "
                     << " mult1:" << (int)((float)timer_ij.elapsed() * 1000) << " Pot:"
                     << (int)((float)timer_p.elapsed() * 1000) << " mult2:" << (int)((float)timer_kij.elapsed() * 1000)
                     << " " << (int)((float)timer_jji.elapsed() * 1000) << " Nnodes: " << N_i << " " << N_j << " "
                     << N_ij << " " << N_p << " " << N_kij << " " << N_jji << " norms " << norm_ij << " " << norm_p
                     << " " << norm_kij << "  " << norm_jji);
}

} // namespace mrchem
