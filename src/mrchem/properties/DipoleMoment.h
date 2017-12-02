#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "TelePrompter.h"

class DipoleMoment {
public:
    DipoleMoment() {
        this->nuclear = Eigen::VectorXd::Zero(3);
        this->electronic = Eigen::VectorXd::Zero(3);
    }
    virtual ~DipoleMoment() { }

    Eigen::VectorXd get() const { return this->nuclear + this->electronic; }
    Eigen::VectorXd& getNuclear() { return this->nuclear; }
    Eigen::VectorXd& getElectronic() { return this->electronic; }

    friend std::ostream& operator<<(std::ostream &o, const DipoleMoment &dipole) {
        Eigen::VectorXd mu = dipole.get();

        int oldPrec = TelePrompter::setPrecision(10);
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
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    Eigen::VectorXd nuclear;
    Eigen::VectorXd electronic;
};

