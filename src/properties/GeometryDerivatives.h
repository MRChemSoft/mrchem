#pragma once

#include "mrchem.h"

namespace mrchem {

class GeometryDerivatives final {
public:
    GeometryDerivatives(int k) {
        this->nuclear = DoubleMatrix::Zero(k, 3);
        this->electronic = DoubleMatrix::Zero(k, 3);
    }
    ~GeometryDerivatives() { }

    DoubleMatrix get() const { return this->nuclear + this->electronic; }
    DoubleMatrix &getNuclear() { return this->nuclear; }
    DoubleMatrix &getElectronic() { return this->electronic; }

    friend std::ostream& operator<<(std::ostream &o, const GeometryDerivatives &geomderiv) {
        DoubleMatrix gd = geomderiv.get();

        int oldPrec = mrcpp::Printer::setPrecision(10);
        double au = gd.norm();
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                   Geometry derivatives                     "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Length of vector:   (au)    " << std::setw(30) << au        <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        for (int k = 0; k < gd.rows(); k++) {
            o<<std::setw(19)<<gd(k,0)<<std::setw(20)<<gd(k,1)<<std::setw(20)<<gd(k,2)<<std::endl;
        }
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    DoubleMatrix nuclear;
    DoubleMatrix electronic;
};

} //namespace mrchem
