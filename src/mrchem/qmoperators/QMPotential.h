#pragma once

#include "QMOperator.h"
#include "QMFunction.h"

/** 
 *  \class QMPotential
 *  \brief Operator defining an external potential
 *
 *  Bla bla bla
 *
 *  \author Stig Rune Jensen
 *  \date 2015
 *  
 */
class QMPotential : public QMFunction, public QMOperator {
public:
    QMPotential(int ab = 1);
    virtual ~QMPotential();

    virtual Orbital* operator() (Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    int adap_build;
    void calcRealPart(Orbital &Vphi, Orbital &phi);
    void calcImagPart(Orbital &Vphi, Orbital &phi, bool adjoint);
};

