#include "LinearResponseSolver.h"
#include "HelmholtzOperatorSet.h"
#include "FockOperator.h"
#include "XCOperator.h"
#include "QMOperatorExp.h"
#include "OrbitalVector.h"
#include "Accelerator.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

LinearResponseSolver::LinearResponseSolver(HelmholtzOperatorSet &helm,
                                           Accelerator *k_x,
                                           Accelerator *k_y)
        : SCF(helm),
          kain_x(k_x),
          kain_y(k_y) {
    this->dynamic = false;
    this->frequency = 0.0;
    this->fOper_0 = 0;
    this->fOper_1 = 0;
    this->fMat_0 = 0;
    this->fMat_x = 0;
    this->fMat_y = 0;
    this->orbitals_0 = 0;
    this->orbitals_x = 0;
    this->orbitals_y = 0;
    this->dOrbitals_x = 0;
    this->dOrbitals_y = 0;
}

LinearResponseSolver::~LinearResponseSolver() {
    this->fOper_0 = 0;
    this->fOper_1 = 0;
    this->fMat_0 = 0;
    this->fMat_x = 0;
    this->fMat_y = 0;
    this->orbitals_0 = 0;
    this->orbitals_x = 0;
    this->orbitals_y = 0;
    this->dOrbitals_x = 0;
    this->dOrbitals_y = 0;
    this->kain_x = 0;
    this->kain_y = 0;
}

void LinearResponseSolver::setupUnperturbed(double prec,
                                            FockOperator &fock,
                                            OrbitalVector &phi,
                                            MatrixXd &F) {
    this->orbitals_0 = &phi;
    this->fMat_0 = &F;
    this->fOper_0 = &fock;
    this->fOper_0->setup(prec);
}

void LinearResponseSolver::clearUnperturbed() {
    this->fOper_0->clear();
    this->fOper_0 = 0;
    this->fMat_0 = 0;
    this->orbitals_0 = 0;
}

void LinearResponseSolver::setup(FockOperator &fock, OrbitalVector &phi_x) {
    if (this->orbitals_0 == 0) MSG_ERROR("Unperturbed system not set up");

    this->dynamic = false;
    this->frequency = 0.0;

    this->fOper_1 = &fock;
    //this->fOper_1->setupUnperturbed();

    this->orbitals_x = &phi_x;
    this->orbitals_y = 0;
    this->dOrbitals_x = new OrbitalVector(phi_x);
    this->dOrbitals_y = 0;

    this->fMat_x = new MatrixXd;
    *this->fMat_x = *this->fMat_0;
    this->fMat_y = 0;
}

void LinearResponseSolver::setup(double omega,
                                 FockOperator &f_oper,
                                 OrbitalVector &phi_x,
                                 OrbitalVector &phi_y) {
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

    if (this->dOrbitals_x != 0) delete this->dOrbitals_x;
    if (this->dOrbitals_y != 0) delete this->dOrbitals_y;

    if (this->fMat_x != 0) delete this->fMat_x;
    if (this->fMat_y != 0) delete this->fMat_y;
    //this->fOper_1->clearUnperturbed();

    this->orbitals_x = 0;
    this->orbitals_y = 0;
    this->dOrbitals_x = 0;
    this->dOrbitals_y = 0;
    this->fOper_1 = 0;

    if (this->kain_x != 0) this->kain_x->clear();
    if (this->kain_y != 0) this->kain_y->clear();

    this->orbError.clear();
    this->property.clear();

    resetPrecision();
}

bool LinearResponseSolver::optimize() {
    MatrixXd &F_0 = *this->fMat_0;
    MatrixXd *F_x = this->fMat_x;
    MatrixXd *F_y = this->fMat_y;
    FockOperator &fock_0 = *this->fOper_0;
    FockOperator &fock_1 = *this->fOper_1;
    OrbitalVector &phi_0 = *this->orbitals_0;
    OrbitalVector *phi_x = this->orbitals_x;
    OrbitalVector *phi_y = this->orbitals_y;
    OrbitalVector *dPhi_x = this->dOrbitals_x;
    OrbitalVector *dPhi_y = this->dOrbitals_y;

    double err_t = 1.0;
    double err_o = 1.0;
    double err_p = 1.0;

    int nIter = 0;
    bool converged = false;
    while(nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycle(nIter);
        adjustPrecision(err_o);

        double orb_prec = getOrbitalPrecision();
        double orb_thrs = getOrbitalThreshold();
        double prop_thrs = getPropertyThreshold();

	fock_1.setup(orb_prec);
	double prop = calcProperty();
	this->property.push_back(prop);

	calcHelmholtzUpdates(phi_x, dPhi_x, fMat_x, false);
	calcHelmholtzUpdates(phi_y, dPhi_y, fMat_y, true);
	this->fOper_1->clear();

        if (nIter > 2) {
            if (this->kain_x != 0) this->kain_x->accelerate(orb_prec, *phi_x, *dPhi_x);
            if (this->kain_y != 0) this->kain_y->accelerate(orb_prec, *phi_y, *dPhi_y);
        }

        err_p = calcPropertyError();
        err_o = calcOrbitalError();
	err_t = calcTotalError();
	this->orbError.push_back(err_t);

	if (phi_x != 0) this->add.inPlace(*phi_x, 1.0, *dPhi_x);
	if (phi_y != 0) this->add.inPlace(*phi_y, 1.0, *dPhi_y);

	if (phi_x != 0) phi_x->setErrors(dPhi_x->getNorms());
	if (phi_y != 0) phi_y->setErrors(dPhi_y->getNorms());

        if (dPhi_x != 0) dPhi_x->clear();
        if (dPhi_y != 0) dPhi_y->clear();

	if (phi_x != 0) printOrbitals(phi_x->getNorms(), *phi_x, 1);
	if (phi_y != 0) printOrbitals(phi_y->getNorms(), *phi_y, 1);

        // Finalize SCF cycle
        timer.stop();
        printProperty();
        printTimer(timer.getWallTime());

        if (err_o < orb_thrs and err_p < prop_thrs) {
            converged = true;
            break;
        }
    }
    if (this->kain_x != 0) this->kain_x->clear();
    if (this->kain_y != 0) this->kain_y->clear();
    printConvergence(converged);
    return converged;
}

void LinearResponseSolver::calcHelmholtzUpdates(OrbitalVector *phi_n,
                                                OrbitalVector *dPhi_n, 
                                                MatrixXd *F, 
                                                bool adjoint) {
    NOT_IMPLEMENTED_ABORT;
    /*
    if (phi_n == 0) return;
    if (dPhi_n == 0) MSG_ERROR("Orbital updates not set up");
    if (F == 0) MSG_ERROR("Fock matrix not set up");
    if (this->orbitals_0 == 0) MSG_ERROR("Unperturbed system not set up");

    this->helmholtz->initialize(this->fMat_0->diagonal());
    OrbitalVector *phi_np1 = new OrbitalVector(*phi_n);
    applyHelmholtzOperators(*phi_np1, *F, *phi_n, adjoint);
    this->add.orthogonalize(*phi_np1, *this->orbitals_0);
    this->add(*dPhi_n, 1.0, *phi_np1, -1.0, *phi_n, true);
    delete phi_np1;
    */
}

OrbitalVector* LinearResponseSolver::setupHelmholtzArguments(FockOperator &fock,
                                                             const Eigen::MatrixXd &M,
                                                             OrbitalVector &phi,
                                                             bool adjoint,
							     bool clearFock) {
    NOT_IMPLEMENTED_ABORT;
}

/*
Orbital* LinearResponseSolver::getHelmholtzArgument(int i, 
                                                    MatrixXd &F,
                                                    OrbitalVector &x,
                                                    bool adjoint) {
    OrbitalVector &phi = *this->orbitals_0;
    Orbital &phi_i = phi.getOrbital(i);
    Orbital &x_i = x.getOrbital(i);
    FockOperator &fock_0 = *this->fOper_0;
    FockOperator &fock_1 = *this->fOper_1;

    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
    MatrixXd LmF = L - F;

    double norm_i = x_i.getSquareNorm();
    if (norm_i > 0.0) {
        norm_i = sqrt(norm_i);
    } else {
        norm_i = 0.0;
    }

    Orbital *part_1 = 0;
    Orbital *part_2 = 0;
    Orbital *part_3 = 0;

    if (norm_i > 0.1*getOrbitalThreshold()) {
        part_1 = fock_0.applyPotential(x_i);
    }

    part_2 = calcMatrixPart(i, LmF, x);

    if (adjoint) {
        part_3 = fock_1.adjoint(phi_i);
    } else {
        part_3 = fock_1(phi_i);
    }

    if (part_1 == 0) part_1 = new Orbital(phi_i);
    if (part_2 == 0) part_2 = new Orbital(phi_i);
    if (part_3 == 0) part_3 = new Orbital(phi_i);

    this->add.orthogonalize(*part_3, phi);

    vector<complex<double> > coefs;
    vector<Orbital *> parts;

    double coef = -1.0/(2.0*pi);
    coefs.push_back(coef);
    coefs.push_back(coef);
    coefs.push_back(coef);
    parts.push_back(part_1);
    parts.push_back(part_2);
    parts.push_back(part_3);

    Orbital *argument = new Orbital(phi_i);
    this->add(*argument, coefs, parts, false);

    if (part_1 != 0) delete part_1;
    if (part_2 != 0) delete part_2;
    if (part_3 != 0) delete part_3;

    return argument;
}
*/

double LinearResponseSolver::calcProperty() {
    FockOperator &fock_1 = *this->fOper_1;
    OrbitalVector &phi_0 = *this->orbitals_0;
    OrbitalVector *phi_x = this->orbitals_x;
    OrbitalVector *phi_y = this->orbitals_y;

    // Static response
    if (phi_y == 0) phi_y = phi_x;

    double property = 0.0;
    QMOperatorExp *h_1 = fock_1.getPerturbationOperator();
    if (h_1 != 0) property += h_1->trace(phi_0, *phi_x, *phi_y);
    return property;
}

double LinearResponseSolver::calcTotalError() const {
    int nOrbs = this->orbitals_0->size();
    VectorXd err_x = VectorXd::Zero(nOrbs);
    VectorXd err_y = VectorXd::Zero(nOrbs);
    if (this->dOrbitals_x != 0) err_x = this->dOrbitals_x->getNorms();
    if (this->dOrbitals_y != 0) err_y = this->dOrbitals_y->getNorms();
    return sqrt(err_x.dot(err_x) + err_y.dot(err_y));
}

double LinearResponseSolver::calcOrbitalError() const {
    double err_o = 0.0;
    if (this->dOrbitals_x != 0) {
        err_o = max(err_o, this->dOrbitals_x->getNorms().maxCoeff());
    }
    if (this->dOrbitals_y != 0) {
        err_o = max(err_o, this->dOrbitals_y->getNorms().maxCoeff());
    }
    return err_o;
}

double LinearResponseSolver::calcPropertyError() const {
    int iter = this->property.size();
    return fabs(getUpdate(this->property, iter, false));
}

void LinearResponseSolver::printProperty() const {
    double prop_0 = 0.0, prop_1 = 0.0;
    int iter = this->property.size();
    if (iter > 1) prop_0 = this->property[iter - 2];
    if (iter > 0) prop_1 = this->property[iter - 1];
    printUpdate(" Property ",  prop_1,  prop_1 -  prop_0);
}

