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
    QMPotential(const QMPotential &pot);
    QMPotential &operator=(const QMPotential &pot);
    virtual ~QMPotential();

    virtual Orbital operator()(Orbital inp);
    virtual Orbital dagger(Orbital inp);

    using QMOperator::operator();
    using QMOperator::dagger;

protected:
    int adap_build;

    mrcpp::FunctionTree<3> *calcRealPart(Orbital &phi, bool dagger);
    mrcpp::FunctionTree<3> *calcImagPart(Orbital &phi, bool dagger);
};

} //namespace mrchem;
