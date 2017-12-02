/* \file NonlinearMaximizer
 *
 *  \date Jan 31, 2013
 *  \author Peter Wind <peter.wind@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif Maximization of nonlinear functional.
 */
#pragma once

#pragma GCC system_header
#include <Eigen/Eigenvalues>
#pragma GCC system_header
#include <Eigen/Dense>

class NonlinearMaximizer {
public:
    NonlinearMaximizer() : N2h(0) { }
    int maximize();

protected:
    int N2h;//size (for orbital localization: N2h = N*(N-1)/2)
    Eigen::MatrixXd hessian;
    Eigen::VectorXd gradient;

    virtual double functional(){ return 0.0; }
    virtual double make_gradient(){ return -1.0; }
    virtual double make_hessian(){ return -1.0; }
    virtual void do_step(Eigen::VectorXd step){ }
};

