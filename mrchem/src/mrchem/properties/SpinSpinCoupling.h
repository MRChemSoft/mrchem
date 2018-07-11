#pragma once

#include "mrchem.h"

namespace mrchem {

class Nucleus;

class SpinSpinCoupling final {
public:
    SpinSpinCoupling(const Nucleus &n_k, const Nucleus &n_l)
        : nuc_K(n_k), nuc_L(n_l) {
        this->diamagnetic = DoubleMatrix::Zero(3,3);
        this->paramagnetic = DoubleMatrix::Zero(3,3);
    }
    ~SpinSpinCoupling() { }

    const Nucleus &getNucleusK() const { return this->nuc_K; }
    const Nucleus &getNucleusL() const { return this->nuc_L; }

    DoubleMatrix get() const { return this->diamagnetic + this->paramagnetic; }
    DoubleMatrix &getDiamagnetic() { return this->diamagnetic; }
    DoubleMatrix &getParamagnetic() { return this->paramagnetic; }

    friend std::ostream& operator<<(std::ostream &o, const SpinSpinCoupling &sscc) {
        DoubleMatrix dia = sscc.diamagnetic;
        DoubleMatrix para = sscc.paramagnetic;
        DoubleMatrix tot = dia + para;

        double isoDShz = dia.trace()/3.0;
        double isoPShz = para.trace()/3.0;
        double isoTShz = isoDShz + isoPShz;

        int oldPrec = mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                    Spin-Spin Coupling Constant             "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(5);
        o<<std::setw(3)  << sscc.getNucleusK().getElement().getSymbol();
        o<<std::setw(26) << sscc.getNucleusK().getCoord()[0];
        o<<std::setw(15) << sscc.getNucleusK().getCoord()[1];
        o<<std::setw(15) << sscc.getNucleusK().getCoord()[2];
        o<<std::endl;
        o<<std::setw(3)  << sscc.getNucleusL().getElement().getSymbol();
        o<<std::setw(26) << sscc.getNucleusL().getCoord()[0];
        o<<std::setw(15) << sscc.getNucleusL().getCoord()[1];
        o<<std::setw(15) << sscc.getNucleusL().getCoord()[2];
        o<<std::endl;
        mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (Hz)       " << std::setw(30) << isoTShz   <<std::endl;
        o<<" Diamagnetic      (Hz)       " << std::setw(30) << isoDShz   <<std::endl;
        o<<" Paramagnetic     (Hz)       " << std::setw(30) << isoPShz   <<std::endl;
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
    const Nucleus nuc_K;
    const Nucleus nuc_L;
    DoubleMatrix diamagnetic;
    DoubleMatrix paramagnetic;
};

} //namespace mrchem
