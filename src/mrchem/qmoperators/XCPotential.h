#pragma once

#include "MWMultiplier.h" 
#include "MWDerivative.h" 
#include "XCOperator.h"

/** 
 *  \class XCPotential
 *  \brief Compute XC potential
 *
 *  TO BE COMPLETED
 *
 *  \author Stig Rune Jensen
 *  \date 2015
 *  
 */
class XCPotential {
public:
    XCPotential(int dO = 1, int dI = 0);
    virtual ~XCPotential() { } 

    int getDerIndex() { return derivativeIndex; };
    int getDerOrder() { return derivativeOrder; };
    void calcPotential(XCFunctional * func,
                       FunctionTree<3> ** xcOutput,
                       Density & density,
                       Density * gradient,
                       DerivativeOperator<3> *derivative,
                       int maxScale);
    
protected:
    
    void calcPotentialLDA(FunctionTree<3> ** xcOutput,
                          Density & density,
                          Density * gradient);
    
    void calcPotentialGGA(bool spin,
                          FunctionTree<3> ** xcOutput,
                          Density & density,
                          Density * gradient,
                          DerivativeOperator<3> *derivative,
                          int maxScale);
    
    FunctionTree<3>* calcDivergence(FunctionTreeVector<3> &inp,
                                    DerivativeOperator<3> *derivative,
                                    int maxScale);
    
    FunctionTree<3>* calcGradDotPotDensVec(FunctionTree<3> &V,
                                           FunctionTreeVector<3> &rho,
                                           DerivativeOperator<3> *derivative,
                                           int maxScale);
    
    FunctionTree<3>* calcPotentialGGA(FunctionTreeVector<3> &xc_funcs,
                                      FunctionTreeVector<3> &dRho_a,
                                      FunctionTreeVector<3> &dRho_b,
                                      DerivativeOperator<3> *derivative,
                                      int maxScale);
                                      
    FunctionTree<3>* potentialFunction;
    int derivativeOrder; // 1=potential, 2=hessian, ....
    int derivativeIndex; // from 0 to k. (0 is pure alpha/total and k
                         // is pure beta/spin). Example:
                         // derivativeOrder = 3 and derivativeIndex=2
                         // can either be the alpha/beta/beta
                         // derivative or the N/S/S derivative)
};

