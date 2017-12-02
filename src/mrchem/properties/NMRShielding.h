#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "Nucleus.h"

class NMRShielding {
public:
    NMRShielding(const Nucleus &n) : nuc(n) {
        this->diamagnetic = Eigen::MatrixXd::Zero(3,3);
        this->paramagnetic = Eigen::MatrixXd::Zero(3,3);
    }
    virtual ~NMRShielding() { }

    const Nucleus& getNucleus() const { return this->nuc; }

    Eigen::MatrixXd get() const { return this->diamagnetic + this->paramagnetic; }
    Eigen::MatrixXd& getDiamagnetic() { return this->diamagnetic; }
    Eigen::MatrixXd& getParamagnetic() { return this->paramagnetic; }

    friend std::ostream& operator<<(std::ostream &o, const NMRShielding &nmr) {
        Eigen::MatrixXd dia = nmr.diamagnetic;
        Eigen::MatrixXd para = nmr.paramagnetic;
        Eigen::MatrixXd tot = dia + para;

        double isoDSppm = dia.trace()/3.0;
        double isoPSppm = para.trace()/3.0;
        double isoTSppm = isoDSppm + isoPSppm;

        int oldPrec = TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                    NMR shielding tensor                    "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(5);
        o<<std::setw(3)  << nmr.getNucleus().getElement().getSymbol();
        o<<std::setw(26) << nmr.getNucleus().getCoord()[0];
        o<<std::setw(15) << nmr.getNucleus().getCoord()[1];
        o<<std::setw(15) << nmr.getNucleus().getCoord()[2];
        o<<std::endl;
        TelePrompter::setPrecision(10);
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
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    const Nucleus nuc;
    Eigen::MatrixXd diamagnetic;
    Eigen::MatrixXd paramagnetic;
};

