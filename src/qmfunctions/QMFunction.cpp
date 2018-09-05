#include "MRCPP/Printer"

#include "QMFunction.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

QMFunction& QMFunction::operator=(const QMFunction &func) {
    if (this != &func) {
        this->re = func.re;
        this->im = func.im;
    }
    return *this;
}

void QMFunction::alloc(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) MSG_FATAL("Function not empty");
        this->re = new mrcpp::FunctionTree<3>(*MRA);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) MSG_FATAL("Function not empty");
        this->im = new mrcpp::FunctionTree<3>(*MRA);
    }
}

void QMFunction::clear(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        this->re = nullptr;
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        this->im = nullptr;
    }
}

void QMFunction::free(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) delete this->re;
        this->re = nullptr;
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) delete this->im;
        this->im = nullptr;
    }
}

int QMFunction::getNNodes(int type) const {
    int nNodes = 0;
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) nNodes += this->real().getNNodes();
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) nNodes += this->imag().getNNodes();
    }
    return nNodes;
}

} //namespace mrchem
