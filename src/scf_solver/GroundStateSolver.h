#pragma once

#include "SCF.h"
#include "SCFEnergy.h"

namespace mrchem {

class GroundStateSolver : public SCF {
public:
    GroundStateSolver(HelmholtzVector &h);
    virtual ~GroundStateSolver();

protected:
    std::vector<SCFEnergy> energy;

    FockOperator *fOper_n;
    Eigen::MatrixXd *fMat_n;

    OrbitalVector *orbitals_n;
    OrbitalVector *orbitals_np1;
    OrbitalVector *dOrbitals_n;

    OrbitalVector *setupHelmholtzArguments(FockOperator &fock,
                                           const Eigen::MatrixXd &M,
                                           OrbitalVector &phi,
                                           bool adjoint = false,
                                           bool clearFock = false);
    void printProperty() const;
    double calcProperty();
    double calcPropertyError() const;

    void localize(FockOperator &fock, Eigen::MatrixXd &F, OrbitalVector &phi);
    void diagonalize(FockOperator &fock, Eigen::MatrixXd &F, OrbitalVector &phi);
    void orthonormalize(FockOperator &fock, Eigen::MatrixXd &F, OrbitalVector &phi);
    Eigen::MatrixXd calcOrthonormalizationMatrix(OrbitalVector &phi);
    Eigen::MatrixXd calcOrthonormalizationMatrix_P(OrbitalVector &phi);
};

/** subclass which defines the particular Gradient and Hessian
 * and other specific functions for a maximization of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 *
 */

#include "NonlinearMaximizer.h"

class RR : public NonlinearMaximizer {
public:
    RR(double prec, OrbitalVector &phi);//make the matrices <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
    const Eigen::MatrixXd &getTotalU() const { return this->total_U; }
protected:
    int N;//number of orbitals
    Eigen::MatrixXd r_i_orig;//<i|R_x|j>,<i|R_y|j>,<i|R_z|j>
    Eigen::MatrixXd r_i ;// rotated  r_i_orig
    Eigen::MatrixXd total_U;// the rotation matrix of the orbitals

    //NB:total_U is not Unitary if the basis set is not orthonormal
    double functional(); //the functional to maximize
    double make_gradient();
    double make_hessian();
    void do_step(Eigen::VectorXd step);
};

} //namespace mrchem
