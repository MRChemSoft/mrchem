#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "FunctionTree.h"
#include "FunctionTreeVector.h"
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
    XCFunctional(bool s, bool e, double thrs = 0.0);
    virtual ~XCFunctional();

    void setDensityCutoff(double thrs) { this->cutoff = thrs; }
    void setFunctional(const std::string &name, double coef = 1.0);

    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength() const { return xc_output_length(this->functional); }

    bool isLDA() const { return (!(this->isGGA() || this->isMetaGGA())); }
    bool isGGA() const { return (xc_is_gga(this->functional)); }
    bool isMetaGGA() const { return (xc_is_metagga(this->functional)); }
    
    bool isSpinSeparated() const { return this->spin; }
    bool needsGamma() const { return (expDerivatives == 0);};

    void evaluate(int k, Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;
    void evalSetup(const int order);

    FunctionTree<3> * calcPotentialGGA(FunctionTree<3> & df_drho, FunctionTree<3> & df_dgamma,
                                       FunctionTreeVector<3> grad_rho, DerivativeOperator<3> *derivative,
                                       int maxScale);
    FunctionTree<3> * calcPotentialGGA(FunctionTree<3> & df_drhoa, FunctionTree<3> & df_dgaa,
                                       FunctionTree<3> & df_dgab, FunctionTreeVector<3> grad_rhoa,
                                       FunctionTreeVector<3> grad_rhob, DerivativeOperator<3> *derivative,
                                       int maxScale);
    FunctionTree<3> * calcPotentialGGA(FunctionTree<3> & df_drho, FunctionTreeVector<3> & df_dgr,
                                       DerivativeOperator<3> *derivative, int maxScale);
 protected:
    FunctionTree<3> * addPotentialContributions(FunctionTreeVector<3> & contributions,
                                                int maxScale);
    FunctionTree<3>* calcDivergence(FunctionTreeVector<3> &inp,
                                    DerivativeOperator<3> *derivative,
                                    int maxScale);
    FunctionTree<3>* calcGradDotPotDensVec(FunctionTree<3> &V,
                                           FunctionTreeVector<3> &rho,
                                           DerivativeOperator<3> *derivative,
                                           int maxScale);

private:
    bool spin;                  ///< Spin polarization
    unsigned int expDerivatives;///< whether gamma-type or explicit derivatives are used
    double cutoff;              ///< Below the cutoff value, the density will be considered zero
    xc_functional functional;   ///< The functional in the XCFun library (struct from xcfun library)

};

