#pragma once

#include <Eigen/Core>

#include "mrchem.h"

/* @class NonlinearMaximizer
 *
 * @date Jan 31, 2013
 * @author Peter Wind <peter.wind@uit.no> \n
 *         CTCC, University of Tromsø
 *
 * @brief Maximization of nonlinear functional
 */

namespace mrchem {

class NonlinearMaximizer {
public:
    NonlinearMaximizer() : N2h(0) { }
    int maximize();

protected:
    int N2h;//size (for orbital localization: N2h = N*(N-1)/2)
    DoubleMatrix hessian;
    DoubleVector gradient;

    virtual double functional() const { return 0.0; }
    virtual double make_gradient() { return -1.0; }
    virtual double make_hessian() { return -1.0; }
    virtual double get_hessian(int i, int j) { return -1.0; }
    virtual void multiply_hessian(DoubleVector &vec, DoubleVector &Hv){}
    virtual void do_step(const DoubleVector &step) { }
};

} //namespace mrchem
