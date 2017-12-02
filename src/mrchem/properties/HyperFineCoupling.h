#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "Nucleus.h"
#include "constants.h"

class HyperFineCoupling {
public:
    HyperFineCoupling(const Nucleus &n) : nuc(n) {
        this->spin_term = Eigen::MatrixXd::Zero(1,1);
        this->fc_term = Eigen::MatrixXd::Zero(1,1);
    }
    virtual ~HyperFineCoupling() { }

    const Nucleus& getNucleus() const { return this->nuc; }

    Eigen::MatrixXd get() const { return this->spin_term + this->fc_term; }
    Eigen::MatrixXd& getSpinTerm() { return this->spin_term; }
    Eigen::MatrixXd& getFermiContactTerm() { return this->fc_term; }

    friend std::ostream& operator<<(std::ostream &o, const HyperFineCoupling &hfc) {
        double fc_term = hfc.fc_term(0,0);
        double spin_term = 1.0/hfc.spin_term(0,0);

        double beta_e = 0.5;        // Bohr magneton
        double beta_N = 0.0002723;  // nuclear magneton
        double g_e = 2.00231931;    // free-electron g-value
        double g_N = hfc.getNucleus().getElement().getGValue();

        double hfcc_g = 0.0;
        double hfcc_hz = 0.0;

        int oldPrec = TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                    HyperFine Coupling Constant             "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(5);
        o<<std::setw(3)  << hfc.getNucleus().getElement().getSymbol();
        o<<std::setw(26) << hfc.getNucleus().getCoord()[0];
        o<<std::setw(15) << hfc.getNucleus().getCoord()[1];
        o<<std::setw(15) << hfc.getNucleus().getCoord()[2];
        o<<std::endl;
        TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" A                    (gauss)" << std::setw(30) << hfcc_g    <<std::endl;
        o<<"                      (MHz)  " << std::setw(30) << hfcc_hz   <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    const Nucleus nuc;
    Eigen::MatrixXd fc_term;
    Eigen::MatrixXd spin_term;
};

