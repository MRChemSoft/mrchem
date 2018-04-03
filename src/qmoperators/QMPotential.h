#pragma once

#include "QMFunction.h"
#include "QMOperator.h"

/** @class QMPotential
 *
 * @brief Operator defining a multiplicative potential
 *
 * Inherits the general features of a complex function from QMFunction and
 * implements the multiplication of this function with an Orbital. The actual
 * function representing the operator needs to be implemented in the derived
 * classes, where the *re and *im FunctionTree pointers should be assigned in
 * the setup() function and deallocated in the clear() function.
 *
 */

namespace mrchem {

class QMPotential : public QMFunction, public QMOperator {
public:
    QMPotential(int adap);
    virtual ~QMPotential();

protected:
    int adap_build;

    virtual Orbital apply(Orbital inp);
    virtual Orbital dagger(Orbital inp);

    mrcpp::FunctionTree<3> *calcRealPart(Orbital &phi, bool dagger);
    mrcpp::FunctionTree<3> *calcImagPart(Orbital &phi, bool dagger);
};

} //namespace mrchem;
