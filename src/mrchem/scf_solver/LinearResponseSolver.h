#pragma once

#include "SCF.h"

class LinearResponseSolver : public SCF {
public:
    LinearResponseSolver(HelmholtzOperatorSet &h,
                         Accelerator *k_x = 0,
                         Accelerator *k_y = 0);
    virtual ~LinearResponseSolver();

    void setupUnperturbed(double prec,
                          FockOperator &fock,
                          OrbitalVector &phi,
                          Eigen::MatrixXd &F); 
    void clearUnperturbed();

    void setup(FockOperator &fock, OrbitalVector &phi_x);
    void setup(double omega,
               FockOperator &fock,
               OrbitalVector &phi_x,
               OrbitalVector &phi_y);
    void clear();

    bool optimize();

protected:
    bool dynamic;
    double frequency;

    FockOperator *fOper_0;
    FockOperator *fOper_1;

    Eigen::MatrixXd *fMat_0;
    Eigen::MatrixXd *fMat_x;
    Eigen::MatrixXd *fMat_y;

    OrbitalVector *orbitals_0;
    OrbitalVector *orbitals_x;
    OrbitalVector *orbitals_y;
    OrbitalVector *dOrbitals_x;
    OrbitalVector *dOrbitals_y;

    Accelerator *kain_x;
    Accelerator *kain_y;

    void calcHelmholtzUpdates(OrbitalVector *phi_n,
                              OrbitalVector *dPhi_n, 
                              Eigen::MatrixXd *F, 
                              bool adjoint);

    OrbitalVector *setupHelmholtzArguments(FockOperator &fock,
                                           const Eigen::MatrixXd &M,
                                           OrbitalVector &phi,
                                           bool adjoint = false,
                                           bool clearFock = false);

    void printProperty() const;

    double calcProperty();
    double calcTotalError() const;
    double calcOrbitalError() const;
    double calcPropertyError() const;
};

