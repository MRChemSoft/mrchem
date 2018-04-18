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

    ComplexMatrix *fMat_n;
    FockOperator  *fOper_n;
    OrbitalVector *orbitals_n;

    OrbitalVector setupHelmholtzArguments(FockOperator &fock,
                                          const ComplexMatrix &M,
                                          OrbitalVector &Phi,
                                          bool adjoint = false,
                                          bool clearFock = false);
    void printProperty() const;
    double calcProperty();
    double calcPropertyError() const;
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
    RR(double prec, OrbitalVector &Phi);//make the matrices <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
    const DoubleMatrix &getTotalU() const { return this->total_U; }
protected:
    int N;//number of orbitals
    DoubleMatrix r_i_orig;//<i|R_x|j>,<i|R_y|j>,<i|R_z|j>
    DoubleMatrix r_i ;// rotated  r_i_orig
    DoubleMatrix total_U;// the rotation matrix of the orbitals

    //NB:total_U is not Unitary if the basis set is not orthonormal
    double functional(); //the functional to maximize
    double make_gradient();
    double make_hessian();
    void do_step(DoubleVector step);
};

} //namespace mrchem
