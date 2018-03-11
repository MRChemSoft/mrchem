#include "MRCPP/Printer"

#include "QMOperator.h"
#include "Orbital.h"

namespace mrchem {

void QMOperator::setApplyPrec(double prec) {
    if (this->apply_prec < 0.0) { 
        this->apply_prec = prec;
    } else if (not isSetup(prec)) {
        MSG_ERROR("Clear operator before setup with different prec!");
    }
}

bool QMOperator::isSetup(double prec) {
    double dPrec = std::abs(this->apply_prec - prec);
    double thrs = mrcpp::MachineZero;
    return (dPrec < thrs) ? true : false;
}

} //namespace mrchem
