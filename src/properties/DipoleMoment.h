#pragma once

#include "mrchem.h"

namespace mrchem {

class DipoleMoment final {
public:
    DipoleMoment() {
        this->nuclear = DoubleVector::Zero(3);
        this->electronic = DoubleVector::Zero(3);
    }
    ~DipoleMoment() { }

    DoubleVector get() const { return this->nuclear + this->electronic; }
    DoubleVector &getNuclear() { return this->nuclear; }
    DoubleVector &getElectronic() { return this->electronic; }

    friend std::ostream& operator<<(std::ostream &o, const DipoleMoment &dipole) {
        DoubleVector mu = dipole.get();

        int oldPrec = mrcpp::Printer::setPrecision(10);
        double au = mu.norm();
        double debye = au/0.393430307;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                         Dipole moment                      "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Length of vector:   (au)    " << std::setw(30) << au        <<std::endl;
        o<<"                     (Debye) " << std::setw(30) << debye     <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<mu(0)<<std::setw(20)<<mu(1)<<std::setw(20)<<mu(2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    DoubleVector nuclear;
    DoubleVector electronic;
};

} //namespace mrchem
