#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "Accelerator.h"
#include "HelmholtzVector.h"
#include "LinearResponseSolver.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/two_electron/FockOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param h: Helmholtz operators
 * @param k_x: Iterative accelerator X orbitals
 * @param k_y: Iterative accelerator Y orbitals
 *
 * SCF solver will NOT take ownership of the HelmholtzVector or the Accelerators,
 * so the original objects must be taken care of externally (do not delete until
 * SCF goes out of scope). Fock matrix, FockOperator and OrbitalVector are not
 * initialized at this stage, so the SCF solver needs to be "setup()" before
 * "optimize()".
 */
LinearResponseSolver::LinearResponseSolver(HelmholtzVector &h, Accelerator *k_x, Accelerator *k_y)
        : SCF(h)
        , dynamic(false)
        , frequency(0.0)
        , fOper_0(nullptr)
        , fOper_1(nullptr)
        , fMat_0(nullptr)
        , fMat_x(nullptr)
        , fMat_y(nullptr)
        , orbitals_0(nullptr)
        , orbitals_x(nullptr)
        , orbitals_y(nullptr)
        , kain_x(k_x)
        , kain_y(k_y) {}

/** @brief Prepare the unperturbed parts of the response solver for optimization
 *
 * @param prec: Precision
 * @param fock: Unperturbed Fock operator
 * @param Phi: Unperturbed orbitals to optimize
 * @param F: Unperturbed Fock matrix
 *
 * The unperturbed part of the response solver remains unchanged during the SCF
 * procedure and for all types and directions of the perturbation, so it needs to
 * be setup only once (unlike the perturbed parts which must be updated in each
 * iteration). SCF solver will NOT take ownership of the input, so these objects
 * must be taken care of externally (do not delete until SCF goes out of scope).
 */
// clang-format off
void LinearResponseSolver::setupUnperturbed(double prec,
                                            FockOperator *fock,
                                            OrbitalVector *Phi,
                                            ComplexMatrix *F) {
    // clang-format on
    this->fOper_0 = fock;
    this->orbitals_0 = Phi;
    this->fMat_0 = F;
    this->fOper_0->setup(prec);
}

/** @brief Clear the unperturbed parts of the response solver after optimization
 *
 * Clear pointers that was set during setupUnperturbed. Solver can be re-used after
 * another setupUnperturbed.
 */
void LinearResponseSolver::clearUnperturbed() {
    this->fOper_0->clear();
    this->fOper_0 = nullptr;
    this->orbitals_0 = nullptr;
    this->fMat_0 = nullptr;
}

/** @brief Prepare solver for optimization, static version
 *
 * @param fock: Perturbed Fock operator (V_1 + h_1)
 * @param X: Perturbed orbitals to optimize (epsilon + omega)
 *
 * SCF solver will NOT take ownership of the input, so these objects must be taken
 * care of externally (do not delete until SCF goes out of scope).
 */
void LinearResponseSolver::setup(FockOperator *fock, OrbitalVector *X) {
    if (this->orbitals_0 == nullptr) MSG_ERROR("Unperturbed system not set up");

    this->dynamic = false;
    this->frequency = 0.0;

    this->fOper_1 = fock;

    this->orbitals_x = X;
    this->orbitals_y = nullptr;

    this->fMat_x = new ComplexMatrix;
    this->fMat_y = nullptr;

    *this->fMat_x = *this->fMat_0;
}

/** @brief Prepare solver for optimization, dynamic version
 *
 * @param omega: Frequency of perturbing field
 * @param fock: Perturbed Fock operator (V_1 + h_1)
 * @param X: Perturbed orbitals to optimize (epsilon + omega)
 * @param Y: Perturbed orbitals to optimize (epsilon - omega)
 *
 * SCF solver will NOT take ownership of the input, so these objects must be taken
 * care of externally (do not delete until SCF goes out of scope).
 */
// clang-format off
void LinearResponseSolver::setup(double omega,
                                 FockOperator *fock,
                                 OrbitalVector *X,
                                 OrbitalVector *Y) {
    // clang-format on
    NOT_IMPLEMENTED_ABORT;
}

/** @brief Clear solver after optimization
 *
 * Clear pointers that was set during setup, and reset the precision parameter
 * (only the current precision orbPrec[0], not the boundary values orbPrec[1,2]).
 * Solver can be re-used after another setup.
 */
void LinearResponseSolver::clear() {
    if (this->fMat_x != nullptr) delete this->fMat_x;
    if (this->fMat_y != nullptr) delete this->fMat_y;

    this->orbitals_x = nullptr;
    this->orbitals_y = nullptr;
    this->fOper_1 = nullptr;

    if (this->kain_x != nullptr) this->kain_x->clear();
    if (this->kain_y != nullptr) this->kain_y->clear();

    this->orbError.clear();
    this->property.clear();

    resetPrecision();
}

/** @brief Run orbital optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm iterates
 * the Sternheimer response equations in integral form. Common implementation for
 * static and dynamic response. Main points of the algorithm:
 *
 * Pre SCF: setup Helmholtz operators with unperturbed energies
 *
 *  1) Setup perturbed Fock operator
 *  2) For X and Y orbitals do:
 *     a) Apply Helmholtz operator on all orbitals
 *     b) Project out occupied space (1 - rho_0)
 *     c) Compute updates and errors
 *     d) Compute KAIN updates
 *  4) Compute property
 *  5) Check for convergence
 *
 */
bool LinearResponseSolver::optimize() {
    ComplexMatrix &F = *this->fMat_0;
    ComplexMatrix &F_x = *this->fMat_x;
    ComplexMatrix &F_y = *this->fMat_y;
    FockOperator &fock_1 = *this->fOper_1;
    OrbitalVector &Phi = *this->orbitals_0;
    OrbitalVector *X_n = this->orbitals_x;
    OrbitalVector *Y_n = this->orbitals_y;
    HelmholtzVector &H = *this->helmholtz;

    double orb_prec = this->orbPrec[0];
    double err_o = 1.0;
    double err_t = 1.0;
    double err_p = 1.0;

    // Setup Helmholtz operators (fixed, based on unperturbed system)
    H.setup(orb_prec, F.real().diagonal());
    ComplexMatrix L = H.getLambdaMatrix();

    // Placeholders for orbital errors
    DoubleVector errors_x = DoubleVector::Zero(Phi.size());
    DoubleVector errors_y = DoubleVector::Zero(Phi.size());

    int nIter = 0;
    bool converged = false;
    while (nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycleHeader(nIter);
        orb_prec = adjustPrecision(err_o);

        // Setup perturbed Fock operator
        fock_1.setup(orb_prec);

        // Iterate X orbitals
        if (X_n != nullptr) {
            // Compute argument: psi_i = V_0*x_i - sum_j [L-F]_ij*x_j + (1 - rho_0)phi_i
            OrbitalVector Psi_n = setupHelmholtzArguments(*X_n, L - F_x, false);

            // Apply Helmholtz operators
            OrbitalVector X_np1 = H(Psi_n);
            orbital::free(Psi_n);

            // Orthogonalize: X_np1 = (1 - rho_0)X_np1
            orbital::orthogonalize(X_np1, Phi);

            // Compute update and errors
            OrbitalVector dX_n = orbital::add(1.0, X_np1, -1.0, *X_n);
            errors_x = orbital::get_norms(dX_n);

            // Compute KAIN update:
            if (this->kain_x != nullptr) {
                orbital::free(X_np1);
                this->kain_x->accelerate(orb_prec, *X_n, dX_n);
                X_np1 = orbital::add(1.0, *X_n, 1.0, dX_n);
            }

            // Free orbitals
            orbital::free(*X_n);
            orbital::free(dX_n);

            // Prepare for next iteration
            *X_n = X_np1;
            orbital::set_errors(*X_n, errors_x);
            printOrbitals(orbital::get_norms(*X_n), *X_n, 0);
        }

        // Iterate Y orbitals
        if (Y_n != nullptr) NOT_IMPLEMENTED_ABORT;

        // Compute property
        double prop = calcProperty();
        this->property.push_back(prop);

        // Compute errors
        err_p = std::abs(getUpdate(this->property, nIter, false));
        err_o = std::max(errors_x.maxCoeff(), errors_y.maxCoeff());
        err_t = std::sqrt(errors_x.dot(errors_x) + errors_y.dot(errors_y));
        converged = checkConvergence(err_o, err_p);
        this->orbError.push_back(err_t);

        // Clear perturbed Fock operator
        fock_1.clear();

        // Finalize SCF cycle
        timer.stop();
        printProperty();
        printCycleFooter(timer.getWallTime());

        if (converged) break;
    }

    // Clear KAIN history and Helmholtz operators
    if (this->kain_x != nullptr) this->kain_x->clear();
    if (this->kain_y != nullptr) this->kain_y->clear();
    H.clear();

    printConvergence(converged);
    return converged;
}

/** @brief Computes the Helmholtz argument for the all orbitals
 *
 * @param Phi_1: Perturbed orbitals
 * @param M: Rotation matrix for second term
 * @param adjoint: Use adjoint of Fock operator
 *
 * Argument contains the unperturbed potential operator acting on perturbed
 * orbital i, and the sum of all perturbed orbitals weighted by the unperturbed
 * Fock matrix, and finally the perturbed Fock operator acting on the unperturbed
 * orbitals (projected out occupied space). The effect of using inexact Helmholtz
 * operators are included in Lambda, which is a diagonal matrix with the actual
 * lambda parameters used in the Helmholtz operators (input matrix M is assumed
 * to be L-F).
 *
 * psi_j = \hat{V}^0 phi^1_i
 *       - \sum_j (\Lambda_{ij}-F^0_{ij})phi^1_j
 *       + (1 - rho^0) V^1 phi^0_i
 *
 */
OrbitalVector LinearResponseSolver::setupHelmholtzArguments(OrbitalVector &Phi_1,
                                                            const ComplexMatrix &M,
                                                            bool adjoint) {
    Timer timer_tot, timer_1(false), timer_2(false), timer_3(false);
    Printer::printHeader(0, "Setting up Helmholtz arguments");
    int oldprec = Printer::setPrecision(5);

    OrbitalVector &Phi_0 = *this->orbitals_0;
    RankZeroTensorOperator &V_0 = this->fOper_0->potential();
    RankZeroTensorOperator V_1 = this->fOper_1->potential() + this->fOper_1->perturbation();

    OrbitalVector out = orbital::param_copy(Phi_1);
    for (int i = 0; i < Phi_1.size(); i++) {
        OrbitalVector arg_parts;
        timer_1.start();
        if (Phi_1[i].norm() > 0.1 * this->orbThrs) {
            Orbital part_1 = V_0(Phi_1[i]);
            arg_parts.push_back(part_1);
        }
        timer_1.stop();
        timer_2.start();
        if (M.row(i).norm() > 0.1 * this->orbThrs) {
            ComplexVector c = M.row(i);
            Orbital part_2 = orbital::linear_combination(c, Phi_1);
            arg_parts.push_back(part_2);
        }
        timer_2.stop();
        timer_3.start();
        Orbital part_3;
        if (adjoint) {
            part_3 = V_1.dagger(Phi_0[i]);
        } else {
            part_3 = V_1(Phi_0[i]);
        }
        part_3.orthogonalize(Phi_0);
        timer_3.stop();
        arg_parts.push_back(part_3);

        ComplexVector c = ComplexVector::Constant(arg_parts.size(), -1.0 / (2.0 * MATHCONST::pi));
        out[i] = orbital::linear_combination(c, arg_parts, -1.0);
        orbital::free(arg_parts);
    }

    Printer::printDouble(0, "            V_0 phi_1", timer_1.getWallTime(), 5);
    Printer::printDouble(0, "            F_0 phi_1", timer_2.getWallTime(), 5);
    Printer::printDouble(0, "(1 - rho_0) V_1 phi_0", timer_3.getWallTime(), 5);

    timer_tot.stop();
    Printer::printFooter(0, timer_tot, 2);
    Printer::setPrecision(oldprec);

    return out;
}

/** @brief Computes the expectation value with the perturbing operator(s) */
double LinearResponseSolver::calcProperty() {
    FockOperator &F_1 = *this->fOper_1;
    OrbitalVector &Phi = *this->orbitals_0;
    OrbitalVector *X = this->orbitals_x;
    OrbitalVector *Y = this->orbitals_y;

    // Static response
    if (Y == nullptr) Y = X;
    return F_1.perturbation().trace(Phi, *X, *Y).real();
}

/** @brief Pretty printing of the computed property with update */
void LinearResponseSolver::printProperty() const {
    double prop_0(0.0), prop_1(0.0);
    int iter = this->property.size();
    if (iter > 1) prop_0 = this->property[iter - 2];
    if (iter > 0) prop_1 = this->property[iter - 1];
    Printer::printHeader(0, "                    Value                  Update      Done ");
    printUpdate(" Property   ", prop_1, prop_1 - prop_0);
    Printer::printSeparator(0, '=');
}

} // namespace mrchem
