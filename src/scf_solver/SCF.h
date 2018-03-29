#pragma once

#include <vector>
#include <string>

#include "qmfunctions.h"

namespace mrchem {

class HelmholtzVector;
class Accelerator;
class FockOperator;

class SCF {
public:
    SCF(HelmholtzVector &h);
    virtual ~SCF();

    virtual bool optimize() = 0;

    double getOrbitalPrecision() const  { return this->orbPrec[0]; }
    double getOrbitalThreshold() const  { return this->orbThrs;    }
    double getPropertyThreshold() const { return this->propThrs;   }

    void setThreshold(double orb_thrs, double prop_thrs);
    void setOrbitalPrec(double init, double final);
    void setMaxIterations(int m_iter) { this->maxIter = m_iter; }

    void setRotation(int iter) { this->rotation = iter; }
    void setCanonical(bool can) { this->canonical = can; }

protected:
    int maxIter;
    int rotation;    ///< number of iterations between localization/diagonalization
    bool canonical;  ///< use localized or canonical orbitals
    double orbThrs;  ///< Convergence threshold orbital update norm
    double propThrs; ///< Convergence threshold property
    double orbPrec[3];

    std::vector<double> orbError;
    std::vector<double> property;

    HelmholtzVector *helmholtz;// Pointer to external object, do not delete!

    bool checkConvergence(double err_o, double err_p) const;
    bool needLocalization(int nIter) const;
    bool needDiagonalization(int nIter) const;

    void adjustPrecision(double error);
    void resetPrecision();

    void printUpdate(const std::string &name, double P, double dP) const;
    double getUpdate(const std::vector<double> &vec, int i, bool absPrec) const;

    void printOrbitals(const DoubleVector &epsilon, const OrbitalVector &Phi, int flag) const;
    void printConvergence(bool converged) const;
    void printCycle(int nIter) const;
    void printTimer(double t) const;
    void printMatrix(int level, const DoubleMatrix &M, const char &name, int pr = 5) const;

    virtual OrbitalVector setupHelmholtzArguments(FockOperator &fock,
                                                  const ComplexMatrix &M,
                                                  OrbitalVector &Phi,
                                                  bool adjoint = false,
                                                  bool clearFock = false) = 0;
};

} //namespace mrchem
