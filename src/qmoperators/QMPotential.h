#pragma once

#include "QMFunction.h"
#include "QMOperator.h"

/** 
 *  \class QMPotential
 *  \brief Operator defining a multiplicative potential
 *
 *  \author Stig Rune Jensen
 *  \date 2015
 *  
 */

namespace mrchem {

class QMPotential : public QMFunction, public QMOperator {
public:
    QMPotential(int adap = 1);
    virtual ~QMPotential();

protected:
    int adap_build;

    virtual Orbital apply(Orbital inp);
    virtual Orbital dagger(Orbital inp);

    using QMOperator::apply;    // Necessary in order to pick up base class 
    using QMOperator::dagger;   // definitions for overloaded functions

    mrcpp::FunctionTree<3> *calcRealPart(Orbital &phi, bool dagger);
    mrcpp::FunctionTree<3> *calcImagPart(Orbital &phi, bool dagger);
};

} //namespace mrchem;
