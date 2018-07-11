#pragma once

#include "mrchem.h"

namespace mrchem {

class Nucleus;

class HyperFineCoupling final {
public:
    HyperFineCoupling(const Nucleus &n) : nuc(n) {
        this->spin_term = DoubleMatrix::Zero(1,1);
        this->fc_term = DoubleMatrix::Zero(1,1);
    }
    ~HyperFineCoupling() { }

    const Nucleus& getNucleus() const { return this->nuc; }

    DoubleMatrix get() const { return this->spin_term + this->fc_term; }
    DoubleMatrix &getSpinTerm() { return this->spin_term; }
    DoubleMatrix &getFermiContactTerm() { return this->fc_term; }

    friend std::ostream& operator<<(std::ostream &o, const HyperFineCoupling &hfc) {
        double fc_term = hfc.fc_term(0,0);
        double spin_term = 1.0/hfc.spin_term(0,0);

        double beta_e = PHYSCONST::beta_e;  // Bohr magneton
        double beta_N = PHYSCONST::beta_N;  // Nuclear magneton
        double g_e    = PHYSCONST::g_e;     // Free-electron g-value
        double g_N    = 0.0; //hfc.getNucleus().getElement().getGValue();

        double hfcc_g = 0.0;
        double hfcc_hz = 0.0;

        int oldPrec = mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                    HyperFine Coupling Constant             "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(5);
        o<<std::setw(3)  << hfc.getNucleus().getElement().getSymbol();
        o<<std::setw(26) << hfc.getNucleus().getCoord()[0];
        o<<std::setw(15) << hfc.getNucleus().getCoord()[1];
        o<<std::setw(15) << hfc.getNucleus().getCoord()[2];
        o<<std::endl;
        mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" A                    (gauss)" << std::setw(30) << hfcc_g    <<std::endl;
        o<<"                      (MHz)  " << std::setw(30) << hfcc_hz   <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    const Nucleus nuc;
    DoubleMatrix fc_term;
    DoubleMatrix spin_term;
};

} //namespace mrchem
