#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "xcfun.h"

/** 
 *  \class XCFunctional
 *  \brief Compute XC functional with XCFun
 *
 *  Interface class for the XCFun library
 *
 *  \author Stig Rune Jensen
 *  \date 2015
 *  
 */
class XCFunctional {
public:
    XCFunctional(bool s, double thrs = 0.0);
    virtual ~XCFunctional();

    void setDensityCutoff(double thrs) { this->cutoff = thrs; }
    void setFunctional(const std::string &name, double coef = 1.0);

    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength() const { return xc_output_length(this->functional); }

    bool isLDA() const { return (!(this->isGGA() && this->isMetaGGA())); }
    bool isGGA() const { return (xc_is_gga(this->functional)); }
    bool isMetaGGA() const { return (xc_is_metagga(this->functional)); }
    
    bool isSpinSeparated() const { return this->spin; }

    void evaluate(int k, Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;

private:
    bool spin;                 ///< Spin polarization
    double cutoff;             ///< Below the cutoff value, the density will be considered zero
    xc_functional functional;  ///< The functional in the XCFun library (struct from xcfun library)

};

