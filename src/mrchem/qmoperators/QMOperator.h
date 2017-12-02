#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "TelePrompter.h"
#include "constants.h"

class Orbital;
class OrbitalVector;

/** \brief Quantum mechanical Hermitian operators

    Base class to handle operators and their application in the QM sense

  */
class QMOperator {
public:
    QMOperator(int ms) : max_scale(ms), apply_prec(-1.0) { }
    QMOperator(const QMOperator &oper) : max_scale(oper.max_scale), apply_prec(oper.apply_prec) { }
    QMOperator& operator=(const QMOperator &inp) { this->apply_prec = inp.apply_prec; return *this; }
    virtual ~QMOperator() { }

    void setApplyPrec(double prec);
    void clearApplyPrec() { this->apply_prec = -1.0; }

    double getApplyPrec() const { return this->apply_prec; }
    int getMaxScale() const { return this->max_scale; }

    virtual void setup(double prec) = 0;
    virtual void clear() = 0;

    virtual Orbital* operator() (Orbital &phi) = 0;
    virtual Orbital* adjoint(Orbital &phi) = 0;

    virtual double operator() (Orbital &phi_i, Orbital &phi_j);
    virtual double adjoint(Orbital &phi_i, Orbital &phi_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    const int max_scale;
    double apply_prec;
};


