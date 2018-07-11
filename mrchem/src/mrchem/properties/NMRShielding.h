#pragma once

#include "mrchem.h"

namespace mrchem {

class Nucleus;

class NMRShielding final {
public:
    NMRShielding(const Nucleus &n) : nuc(n) {
        this->diamagnetic = DoubleMatrix::Zero(3,3);
        this->paramagnetic = DoubleMatrix::Zero(3,3);
    }
    ~NMRShielding() { }

    const Nucleus& getNucleus() const { return this->nuc; }

    DoubleMatrix get() const { return this->diamagnetic + this->paramagnetic; }
    DoubleMatrix &getDiamagnetic() { return this->diamagnetic; }
    DoubleMatrix &getParamagnetic() { return this->paramagnetic; }

    friend std::ostream& operator<<(std::ostream &o, const NMRShielding &nmr) {
        DoubleMatrix dia = nmr.diamagnetic;
        DoubleMatrix para = nmr.paramagnetic;
        DoubleMatrix tot = dia + para;

        double isoDSppm = dia.trace()/3.0;
        double isoPSppm = para.trace()/3.0;
        double isoTSppm = isoDSppm + isoPSppm;

        int oldPrec = mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                    NMR shielding tensor                    "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(5);
        o<<std::setw(3)  << nmr.getNucleus().getElement().getSymbol();
        o<<std::setw(26) << nmr.getNucleus().getCoord()[0];
        o<<std::setw(15) << nmr.getNucleus().getCoord()[1];
        o<<std::setw(15) << nmr.getNucleus().getCoord()[2];
        o<<std::endl;
        mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (ppm)      " << std::setw(30) << isoTSppm  <<std::endl;
        o<<" Diamagnetic      (ppm)      " << std::setw(30) << isoDSppm  <<std::endl;
        o<<" Paramagnetic     (ppm)      " << std::setw(30) << isoPSppm  <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------------- Total ---------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<tot(0,0)<<std::setw(20)<<tot(1,0)<<std::setw(20)<<tot(2,0)<<std::endl;
        o<<std::setw(19)<<tot(1,0)<<std::setw(20)<<tot(1,1)<<std::setw(20)<<tot(2,1)<<std::endl;
        o<<std::setw(19)<<tot(2,0)<<std::setw(20)<<tot(1,2)<<std::setw(20)<<tot(2,2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Diamagnetic ------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<dia(0,0)<<std::setw(20)<<dia(1,0)<<std::setw(20)<<dia(2,0)<<std::endl;
        o<<std::setw(19)<<dia(1,0)<<std::setw(20)<<dia(1,1)<<std::setw(20)<<dia(2,1)<<std::endl;
        o<<std::setw(19)<<dia(2,0)<<std::setw(20)<<dia(1,2)<<std::setw(20)<<dia(2,2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Paramagnetic -----------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<para(0,0)<<std::setw(20)<<para(1,0)<<std::setw(20)<<para(2,0)<<std::endl;
        o<<std::setw(19)<<para(1,0)<<std::setw(20)<<para(1,1)<<std::setw(20)<<para(2,1)<<std::endl;
        o<<std::setw(19)<<para(2,0)<<std::setw(20)<<para(1,2)<<std::setw(20)<<para(2,2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    const Nucleus nuc;
    DoubleMatrix diamagnetic;
    DoubleMatrix paramagnetic;
};

} //namespace mrchem
