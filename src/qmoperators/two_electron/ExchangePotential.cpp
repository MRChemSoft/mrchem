#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotential.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param[in] P Poisson operator (does not take ownership)
 * @param[in] Phi vector of orbitals which define the exchange operator
 */
ExchangePotential::ExchangePotential(mrcpp::PoissonOperator *P, OrbitalVector *Phi, bool s)
        : screen(s)
        , orbitals(Phi)
        , poisson(P) {
    int nOrbs = this->orbitals->size();
    this->tot_norms = DoubleVector::Zero(nOrbs);
    this->part_norms = DoubleMatrix::Zero(nOrbs, nOrbs);
}

/** @brief Prepare operator for application
 *
 * @param[in] prec reqested precision
 *
 * This will NOT precompute the internal exchange between the orbtials defining
 * the operator, which is done explicitly using setupInternal().
 */
void ExchangePotential::setup(double prec) {
    setApplyPrec(prec);

    int nOrbs = this->orbitals->size();
    if (tot_norms.size() != nOrbs) this->tot_norms = DoubleVector::Zero(nOrbs);
    if (part_norms.rows() != nOrbs) this->part_norms = DoubleMatrix::Zero(nOrbs, nOrbs);
    if (part_norms.cols() != nOrbs) this->part_norms = DoubleMatrix::Zero(nOrbs, nOrbs);
}

/** @brief Clears the Exchange Operator
 *
 *  Clears deletes the precomputed exchange contributions.
 */
void ExchangePotential::clear() {
    this->exchange.clear();
    clearApplyPrec();
}

/** @brief Perform a unitary transformation among the precomputed exchange contributions
 *
 * @param[in] U unitary matrix defining the rotation
 */
void ExchangePotential::rotate(const ComplexMatrix &U) {
    if (this->exchange.size() == 0) return;
    this->exchange = orbital::rotate(U, this->exchange, this->apply_prec);

    // NOTE: The following MPI point is currently NOT implemented!
    //
    // the last parameter, 1, means MPI will send only one orbital at a time
    // (because Exchange orbitals can be large for large molecules).
    // OrbitalAdder add(this->apply_prec, this->max_scale, 1);
    // add.rotate(this->exchange, U);
}

/** @brief Applies operator potential
 *
 *  @param[in] inp input orbital
 *
 * The exchange potential is applied to the given orbital. Checks first if this
 * particular exchange contribution has been precomputed.
 */
Orbital ExchangePotential::apply(Orbital inp) {
    if (this->apply_prec < 0.0) {
        MSG_ERROR("Uninitialized operator");
        return inp.paramCopy();
    }
    int i = testPreComputed(inp);
    if (i < 0) {
        println(1, "On-the-fly exchange");
        return calcExchange(inp);
    } else {
        println(1, "Precomputed exchange");
        Orbital out = this->exchange[i].paramCopy();
        qmfunction::deep_copy(out, inp);
        return out;
    }
}

/** @brief Applies the adjoint of the operator
 *  \param[in] inp input orbital
 *
 * NOT IMPLEMENTED
 */
Orbital ExchangePotential::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief Computes the exchange potential on the fly
 *
 *  \param[in] inp input orbital
 *
 * The exchange potential is computed and applied on the fly to the given orbital.
 */
Orbital ExchangePotential::calcExchange(Orbital phi_p) {
    Timer timer;

    double prec = this->apply_prec;
    OrbitalVector &Phi = *this->orbitals;
    mrcpp::PoissonOperator &P = *this->poisson;

    ComplexVector coef_vec(Phi.size());
    QMFunctionVector func_vec;

    OrbitalIterator iter(Phi);
    while (iter.next()) {
        for (int i = 0; i < iter.get_size(); i++) {
            Orbital &phi_i = iter.orbital(i);

            double spin_fac = getSpinFactor(phi_i, phi_p);
            if (std::abs(spin_fac) < mrcpp::MachineZero) continue;

            // compute phi_ip = phi_i^dag * phi_p
            Orbital phi_ip = phi_p.paramCopy();
            qmfunction::multiply(phi_ip, phi_i.dagger(), phi_p, -1.0);

            // compute V_ip = P[phi_ip]
            Orbital V_ip = phi_p.paramCopy();
            if (phi_ip.function().hasReal()) {
                V_ip.function().alloc(NUMBER::Real);
                mrcpp::apply(prec, V_ip.function().real(), P, phi_ip.function().real());
            }
            if (phi_ip.function().hasImag()) {
                V_ip.function().alloc(NUMBER::Imag);
                mrcpp::apply(prec, V_ip.function().imag(), P, phi_ip.function().imag());
            }
            phi_ip.release();

            // compute phi_iip = phi_i * V_ip
            Orbital phi_iip = phi_p.paramCopy();
            qmfunction::multiply(phi_iip, phi_i, V_ip, -1.0);

            coef_vec(i) = spin_fac / phi_i.squaredNorm();
            func_vec.push_back(phi_iip);
        }
    }

    // compute ex_p = sum_i c_i*phi_iip
    Orbital ex_p = phi_p.paramCopy();
    qmfunction::linear_combination(ex_p, coef_vec, func_vec, -1.0);

    timer.stop();
    double n = ex_p.function().getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied exchange", n, t);

    return ex_p;
}

/** @brief precomputes the exchange potential
 *
 *  @param[in] phi_p input orbital
 *
 * The exchange potential is (pre)computed among the orbitals that define the operator
 */
void ExchangePotential::setupInternal(double prec) {
    setApplyPrec(prec);

    if (this->exchange.size() != 0) MSG_ERROR("Exchange not properly cleared");

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Ex = this->exchange;

    Timer timer;
    // Diagonal must come first because it's NOT in-place
    for (int i = 0; i < Phi.size(); i++) calcInternal(i);

    // Off-diagonal must come last because it IS in-place
    OrbitalIterator iter(Phi, true); //symmetric iterator
    Orbital ex_rcv;
    while (iter.next(1)) { //one orbital at the time
        if (iter.get_size() > 0) {
            Orbital &phi_i = iter.orbital(0);
            int idx = iter.idx(0);
            for (int j = 0; j < Phi.size(); j++) {
                if (mpi::my_orb(phi_i) and j <= idx) continue; //compute only i<j for own block
                Orbital &phi_j = (*this->orbitals)[j];
                if (mpi::my_orb(phi_j)) calcInternal(idx, j, phi_i, phi_j);
            }
            //must send exchange_i to owner and receive exchange computed by other
            if (iter.get_step(0) and not mpi::my_orb(phi_i))
                mpi::send_function(Ex[idx].function(), phi_i.rankID(), idx, mpi::comm_orb);

            if (iter.get_sent_size()) {
                //get exchange from where we sent orbital to
                int idx_sent = iter.get_idx_sent(0);
                int sent_rank = iter.get_rank_sent(0);
                mpi::recv_function(ex_rcv.function(), sent_rank, idx_sent, mpi::comm_orb);
                Ex[idx_sent].add(1.0, ex_rcv);
            }

            if (not iter.get_step(0) and not mpi::my_orb(phi_i))
                mpi::send_function(Ex[idx].function(), phi_i.rankID(), idx, mpi::comm_orb);
            if (not mpi::my_orb(Ex[idx])) Ex[idx].freeFunctions();
        } else {
            if (iter.get_sent_size()) { //must receive exchange computed by other
                //get exchange from where we sent orbital to
                int idx_sent = iter.get_idx_sent(0);
                int sent_rank = iter.get_rank_sent(0);
                mpi::recv_function(ex_rcv.function(), sent_rank, idx_sent, mpi::comm_orb);
                Ex[idx_sent].add(1.0, ex_rcv);
            }
        }
        ex_rcv.freeFunctions();
    }

    int n = 0;
    // Collect info from the calculation
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) this->tot_norms(i) = Ex[i].norm();
        n += Ex[i].function().getNNodes(NUMBER::Total);
    }

    mpi::allreduce_vector(this->tot_norms, mpi::comm_orb);  //to be checked
    mpi::allreduce_matrix(this->part_norms, mpi::comm_orb); //to be checked

    timer.stop();
    double t = timer.getWallTime();
    Printer::printTree(0, "Hartree-Fock exchange", n, t);
}

/** @brief Computes the diagonal part of the internal exchange potential
 *
 *  \param[in] i orbital index
 *
 * The diagonal term K_ii is computed.
 */
void ExchangePotential::calcInternal(int i) {
    Orbital &phi_i = (*this->orbitals)[i];

    if (mpi::my_orb(phi_i)) {
        double prec = std::min(getScaledPrecision(i, i), 1.0e-1);
        mrcpp::PoissonOperator &P = *this->poisson;

        // compute phi_ii = phi_i^dag * phi_i
        Orbital phi_ii = phi_i.paramCopy();
        qmfunction::multiply(phi_ii, phi_i.dagger(), phi_i, prec);

        // compute V_ii = P[phi_ii]
        Orbital V_ii = phi_i.paramCopy();
        if (phi_ii.function().hasReal()) {
            V_ii.function().alloc(NUMBER::Real);
            mrcpp::apply(prec, V_ii.function().real(), P, phi_ii.function().real());
        }
        if (phi_ii.function().hasImag()) {
            V_ii.function().alloc(NUMBER::Imag);
            mrcpp::apply(prec, V_ii.function().imag(), P, phi_ii.function().imag());
        }
        phi_ii.release();

        // compute phi_iii = phi_i * V_ii
        Orbital phi_iii = phi_i.paramCopy();
        qmfunction::multiply(phi_iii, phi_i, V_ii, prec);
        phi_iii.rescale(1.0 / phi_i.squaredNorm());
        this->part_norms(i, i) = phi_iii.norm();
        this->exchange.push_back(phi_iii);
    } else {
        //put empty orbital to fill the exchange vector
        Orbital phi_iii = phi_i.paramCopy();
        this->exchange.push_back(phi_iii);
    }
}

/** @brief computes the off-diagonal part of the exchange potential
 *
 *  \param[in] i first orbital index
 *  \param[in] j second orbital index
 *
 * The off-diagonal terms K_ij and K_ji are computed.
 */
void ExchangePotential::calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j) {
    mrcpp::PoissonOperator &P = *this->poisson;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Ex = this->exchange;

    if (i == j) MSG_FATAL("Cannot handle diagonal term");
    if (Ex.size() != Phi.size()) MSG_FATAL("Size mismatch");
    if (phi_i.function().hasImag() or phi_j.function().hasImag()) MSG_FATAL("Orbitals must be real");

    double i_fac = getSpinFactor(phi_i, phi_j);
    double j_fac = getSpinFactor(phi_j, phi_i);

    double thrs = mrcpp::MachineZero;
    if (std::abs(i_fac) < thrs or std::abs(j_fac) < thrs) {
        this->part_norms(i, j) = 0.0;
        return;
    }

    // set correctly scaled precision for components ij and ji
    double prec = std::min(getScaledPrecision(i, j), getScaledPrecision(j, i));
    if (prec > 1.0e00) return;     // orbital does not contribute within the threshold
    prec = std::min(prec, 1.0e-1); // very low precision does not work properly

    // compute phi_ij = phi_i^dag * phi_j (dagger NOT used, orbitals must be real!)
    Orbital phi_ij = phi_i.paramCopy();
    qmfunction::multiply(phi_ij, phi_i, phi_j, prec);

    // compute V_ij = P[phi_ij]
    Orbital V_ij = phi_i.paramCopy();
    if (phi_ij.function().hasReal()) {
        V_ij.function().alloc(NUMBER::Real);
        mrcpp::apply(prec, V_ij.function().real(), P, phi_ij.function().real());
    }
    if (phi_ij.function().hasImag()) {
        MSG_FATAL("Orbitals must be real");
        V_ij.function().alloc(NUMBER::Imag);
        mrcpp::apply(prec, V_ij.function().imag(), P, phi_ij.function().imag());
    }
    phi_ij.release();

    // compute phi_jij = phi_j * V_ij
    Orbital phi_jij = phi_j.paramCopy();
    qmfunction::multiply(phi_jij, phi_j, V_ij, prec);
    phi_jij.rescale(1.0 / phi_j.squaredNorm());
    this->part_norms(j, i) = phi_jij.norm();

    // compute phi_iij = phi_i * V_ij
    Orbital phi_iij = phi_i.paramCopy();
    qmfunction::multiply(phi_iij, phi_i, V_ij, prec);
    phi_iij.rescale(1.0 / phi_i.squaredNorm());
    this->part_norms(i, j) = phi_iij.norm();
    V_ij.release();

    // compute x_i += phi_jij
    Ex[i].add(i_fac, phi_jij);
    phi_jij.release();

    // compute x_j += phi_iij
    Ex[j].add(j_fac, phi_iij);
    phi_iij.release();
}

/** @brief Test if a given contribution has been precomputed
 *
 * @param[in] phi_p orbital for which the check is performed
 *
 * If the given contribution has been precomputed, it is simply copied,
 * without additional recalculation.
 */
int ExchangePotential::testPreComputed(Orbital phi_p) const {
    const OrbitalVector &Phi = *this->orbitals;
    const OrbitalVector &Ex = this->exchange;

    int out = -1;
    if (Ex.size() == Phi.size()) {
        for (int i = 0; i < Phi.size(); i++) {
            if (&Phi[i].function().real() == &phi_p.function().real() and
                &Phi[i].function().imag() == &phi_p.function().imag()) {
                out = i;
                break;
            }
        }
    }
    return out;
}

/** @brief scale the relative precision based on norm
 *
 * The internal norms are saved between SCF iterations so that they can
 * be used to estimate the size of the different contributions to the total
 * exchange. The relative precision of the Poisson terms is scaled to give a
 * consistent _absolute_ pecision in the final output.
 */
double ExchangePotential::getScaledPrecision(int i, int j) const {
    double scaled_prec = this->apply_prec;
    if (this->screen) {
        double tNorm = this->tot_norms(i);
        double pNorm = std::max(this->part_norms(i, j), this->part_norms(j, i));
        if (tNorm > 0.0) scaled_prec *= tNorm / pNorm;
    }
    return scaled_prec;
}

/** @brief determines the exchange factor to be used in the calculation of the exact exchange
  *
  * @param [in] orb input orbital to which K is applied
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
    if (phi_j.spin() == SPIN::Paired)
        out = 1.0;
    else if (phi_i.spin() == SPIN::Paired)
        out = 0.5;
    else if (phi_i.spin() == phi_j.spin())
        out = 1.0;
    return out;
}

} //namespace mrchem
