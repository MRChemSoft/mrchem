#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "LinearResponseSolver.h"
#include "HelmholtzVector.h"
#include "FockOperator.h"
#include "Accelerator.h"
#include "Orbital.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

LinearResponseSolver::LinearResponseSolver(HelmholtzVector &h,
                                           Accelerator *k_x,
                                           Accelerator *k_y)
        : SCF(h),
          dynamic(false),
          frequency(0.0),
          fOper_0(nullptr),
          fOper_1(nullptr),
          fMat_0(nullptr),
          fMat_x(nullptr),
          fMat_y(nullptr),
          orbitals_0(nullptr),
          orbitals_x(nullptr),
          orbitals_y(nullptr),
          kain_x(k_x),
          kain_y(k_y) {
}

void LinearResponseSolver::setupUnperturbed(double prec,
                                            FockOperator *fock,
                                            OrbitalVector *Phi,
                                            ComplexMatrix *F) {
    this->fOper_0 = fock;
    this->orbitals_0 = Phi;
    this->fMat_0 = F;
    this->fOper_0->setup(prec);
}

void LinearResponseSolver::clearUnperturbed() {
    this->fOper_0->clear();
    this->fOper_0 = nullptr;
    this->orbitals_0 = nullptr;
    this->fMat_0 = nullptr;
}

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

void LinearResponseSolver::setup(double omega,
                                 FockOperator *fock,
                                 OrbitalVector *X,
                                 OrbitalVector *Y) {
    NOT_IMPLEMENTED_ABORT;
    /*
    if (this->orbitals_0 == 0) MSG_ERROR("Unperturbed system not set up");

    this->dynamic = true;
    this->frequency = omega;

    this->fOper_1 = &f_oper;
    XCOperator *xc = this->fOper_1->getXCOperator();
    if (xc != 0) xc->setupUnperturbed();

    this->phi_x = &orbs_x;
    this->phi_y = &orbs_y;
    this->dPhi_x = new OrbitalVector("dOrbs_x", orbs_x);
    this->dPhi_y = new OrbitalVector("dOrbs_y", orbs_y);

    bool im = this->fOper_1->getPerturbationOperator(0).isImaginary();
    this->phi_x->setImaginary(im);
    this->phi_y->setImaginary(im);

    int nOrbs = this->phi_0->size();
    MatrixXd Omega = MatrixXd::Identity(nOrbs,nOrbs);
    Omega *= omega;
    this->fMat_x = new MatrixXd;
    this->fMat_y = new MatrixXd;
    *this->fMat_x = *this->fMat_0 + Omega;
    *this->fMat_y = *this->fMat_0 - Omega;
    */
}

void LinearResponseSolver::clear() {
    this->clearUnperturbed();

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

bool LinearResponseSolver::optimize() {
    ComplexMatrix &F = *this->fMat_0;
    ComplexMatrix *F_x = this->fMat_x;
    ComplexMatrix *F_y = this->fMat_y;
    FockOperator &fock_0 = *this->fOper_0;
    FockOperator &fock_1 = *this->fOper_1;
    OrbitalVector &Phi = *this->orbitals_0;
    OrbitalVector *X_n = this->orbitals_x;
    OrbitalVector *Y_n = this->orbitals_y;
    HelmholtzVector &H = *this->helmholtz;

    double orb_prec = this->orbPrec[0];
    double err_o = 1.0;
    double err_t = 1.0;
    double err_p = 1.0;

    int nIter = 0;
    bool converged = false;
    while(nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycleHeader(nIter);
        orb_prec = adjustPrecision(err_o);

        fock_1.setup(orb_prec);
        double prop = calcProperty();
        this->property.push_back(prop);

        DoubleVector errors_x = DoubleVector::Zero(Phi.size());
        DoubleVector errors_y = DoubleVector::Zero(Phi.size());

        H.setup(orb_prec, F.real().diagonal());
        if (X_n != nullptr) {
            ComplexMatrix L = H.getLambdaMatrix();
            OrbitalVector Psi_n = setupHelmholtzArguments(*X_n, L-F, false);
            OrbitalVector X_np1 = H(Psi_n);
            orbital::orthogonalize(X_np1, Phi);
            OrbitalVector dX_n = orbital::add(1.0, X_np1, -1.0, *X_n);
            errors_x = orbital::get_norms(dX_n);
            if (this->kain_x != nullptr) {
                orbital::free(X_np1);
                this->kain_x->accelerate(orb_prec, *X_n, dX_n);
                X_np1 = orbital::add(1.0, *X_n, 1.0, dX_n);
            }
            orbital::free(*X_n);
            orbital::free(dX_n);
            *X_n = X_np1;
            orbital::set_errors(*X_n, errors_x);
            printOrbitals(orbital::get_norms(*X_n), *X_n, 0);
        }
        if (Y_n != nullptr) {
            NOT_IMPLEMENTED_ABORT;
        }
        H.clear();
        fock_1.clear();

        err_p = calcPropertyError();
        err_o = std::max(errors_x.maxCoeff(), errors_y.maxCoeff());
        err_t = std::sqrt(errors_x.dot(errors_x) + errors_y.dot(errors_y));
        this->orbError.push_back(err_t);
        converged = checkConvergence(err_o, err_p);

        // Finalize SCF cycle
        timer.stop();
        printProperty();
        printCycleFooter(timer.getWallTime());

        if (converged) break;
    }

    if (this->kain_x != nullptr) this->kain_x->clear();
    if (this->kain_y != nullptr) this->kain_y->clear();

    printConvergence(converged);
    return converged;
}

OrbitalVector LinearResponseSolver::setupHelmholtzArguments(OrbitalVector &Phi_1,
                                                            const ComplexMatrix &M,
                                                            bool adjoint) {
    Timer timer_tot;
    Printer::printHeader(0, "Setting up Helmholtz arguments");
    int oldprec = Printer::setPrecision(5);

    OrbitalVector &Phi_0 = *this->orbitals_0;
    RankZeroTensorOperator &V_0 = this->fOper_0->potential();
    RankZeroTensorOperator V_1 = this->fOper_1->potential() + this->fOper_1->perturbation();

    OrbitalVector out = orbital::param_copy(Phi_1);
    for (int i = 0; i < Phi_1.size(); i++) {
        OrbitalVector arg_parts;
        if (Phi_1[i].norm() > 0.1*this->orbThrs) {
            Orbital part_1 = V_0(Phi_1[i]);
            arg_parts.push_back(part_1);
        }
        if (M.row(i).norm() > 0.1*this->orbThrs) {
            ComplexVector c = M.row(i);
            Orbital part_2 = orbital::multiply(c, Phi_1);
            arg_parts.push_back(part_2);
        }
        Orbital part_3;
        if (adjoint) {
            part_3 = V_1.dagger(Phi_0[i]);
        } else {
            part_3 = V_1(Phi_0[i]);
        }
        part_3.orthogonalize(Phi_0);
        arg_parts.push_back(part_3);

        ComplexVector c = ComplexVector::Constant(arg_parts.size(), -1.0/(2.0*MATHCONST::pi));
        out[i] = orbital::multiply(c, arg_parts, -1.0);
        orbital::free(arg_parts);
    }

    timer_tot.stop();
    Printer::printFooter(0, timer_tot, 2);
    Printer::setPrecision(oldprec);

    return out;
}

double LinearResponseSolver::calcProperty() {
    FockOperator &F_1 = *this->fOper_1;
    OrbitalVector &Phi = *this->orbitals_0;
    OrbitalVector *X = this->orbitals_x;
    OrbitalVector *Y = this->orbitals_y;

    // Static response
    if (Y == nullptr) Y = X;
    return F_1.perturbation().trace(Phi, *X, *Y).real();
}

double LinearResponseSolver::calcPropertyError() const {
    int iter = this->property.size();
    return std::abs(getUpdate(this->property, iter, false));
}

void LinearResponseSolver::printProperty() const {
    double prop_0(0.0), prop_1(0.0);
    int iter = this->property.size();
    if (iter > 1) prop_0 = this->property[iter - 2];
    if (iter > 0) prop_1 = this->property[iter - 1];
    printUpdate(" Property ",  prop_1,  prop_1 -  prop_0);
}

}  //namespace mrchem
