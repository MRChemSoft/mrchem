#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotentialD2.h"
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
 */
ExchangePotentialD2::ExchangePotentialD2(PoissonOperator_p P,
                                         OrbitalVector_p Phi,
                                         OrbitalVector_p X,
                                         OrbitalVector_p Y,
                                         double prec)
        : ExchangePotential(P, Phi, prec)
        , orbitals_x(X)
        , orbitals_y(Y) {
    if (X == Y) useOnlyX = true;
}

/** @brief Save all orbitals in Bank, so that they can be accessed asynchronously */
void ExchangePotentialD2::setupBank() {
    if (mpi::bank_size < 1) return;

    Timer timer;
    mpi::barrier(mpi::comm_orb);
    OrbitalVector &Phi = *this->orbitals;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) PhiBank.put_orb(i, Phi[i]);
    }
    OrbitalVector &X = *this->orbitals_x;
    for (int i = 0; i < X.size(); i++) {
        if (mpi::my_orb(X[i])) XBank.put_orb(i, X[i]);
    }
    OrbitalVector &Y = *this->orbitals_y;
    for (int i = 0; i < Y.size(); i++) {
        if (mpi::my_orb(Y[i])) YBank.put_orb(i, Y[i]);
    }
    mpi::barrier(mpi::comm_orb);
    mrcpp::print::time(4, "Setting up exchange bank", timer);
}

/** @brief Clears the Exchange Operator
 *
 *  Clears the orbital bank accounts.
 */
void ExchangePotentialD2::clearBank() {
    PhiBank.clear();
    XBank.clear();
    YBank.clear();
}

/** @brief Apply exchange operator to given orbital
 *
 *  @param[in] phi_p input orbital
 *
 * The D2 operator has to be applied on-the-fly, e.i. no
 * pre-computed exchange contributions are available.
 */
Orbital ExchangePotentialD2::apply(Orbital phi_p) {
    if (this->apply_prec < 0.0) {
        MSG_ERROR("Uninitialized operator");
        return phi_p.paramCopy();
    }

    int id_shift = 10000;

    Timer timer;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    double prec = this->apply_prec;
    // use fixed exchange_prec if set explicitly, otherwise use setup prec
    double precf = (this->exchange_prec > 0.0) ? this->exchange_prec : prec;
    // adjust precision since we sum over orbitals
    precf /= std::min(10.0, std::sqrt(1.0 * Phi.size()));

    QMFunctionVector func_vec;
    std::vector<ComplexDouble> coef_vec;
    for (int i = 0; i < Phi.size(); i++) {
        Orbital &phi_i = Phi[i];
        Orbital &x_i = X[i];
        Orbital &y_i = Y[i];

        if (not mpi::my_orb(phi_i)) PhiBank.get_orb(i + id_shift, phi_i, 1);
        if (not mpi::my_orb(x_i)) PhiBank.get_orb(i + 2 * id_shift, x_i, 1);
        if (not mpi::my_orb(y_i)) PhiBank.get_orb(i + 3 * id_shift, y_i, 1);

        double spin_fac = getSpinFactor(phi_i, phi_p);
        if (std::abs(spin_fac) >= mrcpp::MachineZero) {
            Orbital ex_xip = phi_p.paramCopy();
            Orbital ex_iyp = phi_p.paramCopy();
            calcExchange_kij(precf, x_i, phi_i, phi_p, ex_xip);
            calcExchange_kij(precf, phi_i, y_i, phi_p, ex_iyp);
            func_vec.push_back(ex_xip);
            func_vec.push_back(ex_iyp);
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
        }
        if (not mpi::my_orb(phi_i)) phi_i.free(NUMBER::Total);
        if (not mpi::my_orb(x_i)) x_i.free(NUMBER::Total);
        if (not mpi::my_orb(y_i)) y_i.free(NUMBER::Total);
    }

    // compute out_p = sum_i c_i*(ex_xip + ex_iyp)
    Orbital out_p = phi_p.paramCopy();
    Eigen::Map<ComplexVector> coefs(coef_vec.data(), coef_vec.size());
    qmfunction::linear_combination(out_p, coefs, func_vec, prec);
    print_utils::qmfunction(3, "Applied exchange", out_p, timer);
    return out_p;
}

/** @brief Apply adjoint exchange operator to given orbital
 *
 *  @param[in] phi_p input orbital
 *
 * The D2 operator has to be applied on-the-fly, e.i. no
 * pre-computed exchange contributions are available.
 */
Orbital ExchangePotentialD2::dagger(Orbital phi_p) {
    if (this->apply_prec < 0.0) {
        MSG_ERROR("Uninitialized operator");
        return phi_p.paramCopy();
    }

    int id_shift = 10000;

    Timer timer;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    double prec = this->apply_prec;
    // use fixed exchange_prec if set explicitly, otherwise use setup prec
    double precf = (this->exchange_prec > 0.0) ? this->exchange_prec : prec;
    // adjust precision since we sum over orbitals
    precf /= std::min(10.0, std::sqrt(1.0 * Phi.size()));

    QMFunctionVector func_vec;
    std::vector<ComplexDouble> coef_vec;
    for (int i = 0; i < Phi.size(); i++) {
        Orbital &phi_i = Phi[i];
        Orbital &x_i = X[i];
        Orbital &y_i = Y[i];

        if (not mpi::my_orb(phi_i)) PhiBank.get_orb(i + id_shift, phi_i, 1);
        if (not mpi::my_orb(x_i)) PhiBank.get_orb(i + 2 * id_shift, x_i, 1);
        if (not mpi::my_orb(y_i)) PhiBank.get_orb(i + 3 * id_shift, y_i, 1);

        double spin_fac = getSpinFactor(phi_i, phi_p);
        if (std::abs(spin_fac) >= mrcpp::MachineZero) {
            Orbital ex_ixp = phi_p.paramCopy();
            Orbital ex_yip = phi_p.paramCopy();
            calcExchange_kij(precf, phi_i, x_i, phi_p, ex_ixp);
            calcExchange_kij(precf, y_i, phi_i, phi_p, ex_yip);
            func_vec.push_back(ex_ixp);
            func_vec.push_back(ex_yip);
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
            coef_vec.push_back(spin_fac / phi_i.squaredNorm());
        }
        if (not mpi::my_orb(phi_i)) phi_i.free(NUMBER::Total);
        if (not mpi::my_orb(x_i)) x_i.free(NUMBER::Total);
        if (not mpi::my_orb(y_i)) y_i.free(NUMBER::Total);
    }

    // compute ex_p = sum_i c_i*(ex_ixp + ex_yip)
    Orbital ex_p = phi_p.paramCopy();
    Eigen::Map<ComplexVector> coefs(coef_vec.data(), coef_vec.size());
    qmfunction::linear_combination(ex_p, coefs, func_vec, prec);
    print_utils::qmfunction(3, "Applied exchange", ex_p, timer);
    return ex_p;
}

} // namespace mrchem
