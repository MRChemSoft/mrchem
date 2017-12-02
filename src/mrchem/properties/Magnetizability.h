#pragma once

#pragma GCC system_header
#include <Eigen/Core>

class Magnetizability {
public:
    Magnetizability() {
        this->diamagnetic = Eigen::MatrixXd::Zero(3,3);
        this->paramagnetic = Eigen::MatrixXd::Zero(3,3);
    }
    virtual ~Magnetizability() { }

    Eigen::MatrixXd get() const { return this->diamagnetic + this->paramagnetic; }
    Eigen::MatrixXd& getDiamagnetic() { return this->diamagnetic; }
    Eigen::MatrixXd& getParamagnetic() { return this->paramagnetic; }

    friend std::ostream& operator<<(std::ostream &o, const Magnetizability &mag) {
        Eigen::MatrixXd dia = mag.diamagnetic;
        Eigen::MatrixXd para = mag.paramagnetic;
        Eigen::MatrixXd tot = dia + para;

        double w_au = 0.0;  // Only static magnetizability
        double isoDMau = dia.trace()/3.0;
        double isoPMau = para.trace()/3.0;
        double isoTMau = isoDMau + isoPMau;

        double isoDMsi = isoDMau*78.9451185;
        double isoPMsi = isoPMau*78.9451185;
        double isoTMsi = isoTMau*78.9451185;

        int oldPrec = TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                   Magnetizability tensor                   "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Frequency:       (au)       " << std::setw(30) << w_au      <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (au)       " << std::setw(30) << isoTMau   <<std::endl;
        o<<" Diamagnetic      (au)       " << std::setw(30) << isoDMau   <<std::endl;
        o<<" Paramagnetic     (au)       " << std::setw(30) << isoPMau   <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (SI)       " << std::setw(30) << isoTMsi   <<std::endl;
        o<<" Diamagnetic      (SI)       " << std::setw(30) << isoDMsi   <<std::endl;
        o<<" Paramagnetic     (SI)       " << std::setw(30) << isoPMsi   <<std::endl;
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
    Eigen::MatrixXd diamagnetic;
    Eigen::MatrixXd paramagnetic;
};

