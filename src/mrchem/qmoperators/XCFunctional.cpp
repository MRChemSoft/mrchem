#include "XCFunctional.h"
#include "TelePrompter.h"
#include "constants.h"

using namespace std;
using namespace Eigen;

/** \brief creator 
 *
 * Initializes some parameters in the xcfun library
 *
 */
XCFunctional::XCFunctional(bool s, double thrs)
        : spin(s), cutoff(thrs) {
    this->functional = xc_new_functional();
    if (this->spin) {
        xc_set(this->functional, "XC_A_B", 1);
    } else {
        xc_set(this->functional, "XC_N", 1);
    }
}

XCFunctional::~XCFunctional() {
    xc_free_functional(this->functional);
}

void XCFunctional::setFunctional(const string &name, double coef) {
    xc_set(this->functional, name.c_str(), coef);
}

/** \breif Evaluates XC functional and derivatives
 *
 * Computes the alpha and beta exchange-correlation functionals and
 * their derivatives.  The electronic density (total/alpha/beta) are
 * given as input. Results are then stored in the xcfun output
 * functions. Higher order derivatives can be computed changing the parameter k. 
 *
 * XCFunctional output (with k=1):
 *
 * LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho}\right) \f$
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \rho_x},
 *  \frac{\partial F_{xc}}{\partial \rho_y},
 *  \frac{\partial F_{xc}}{\partial \rho_z}\right) \f$
 *
 * Spin LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta}\right) \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\beta}
 *  \right) \f$
 */
void XCFunctional::evaluate(int k, MatrixXd &inp, MatrixXd &out) const {
    if (inp.cols() != getInputLength()) MSG_ERROR("Invalid input");

    int nInp = getInputLength();
    int nOut = getOutputLength();
    int nPts = inp.rows();
    out = MatrixXd::Zero(nPts, nOut);

    double *iDat = new double[nInp];
    double *oDat = new double[nOut];

    for (int i = 0; i < nPts; i++) {
        if (inp(i,0) > this->cutoff) {
            for (int j = 0; j < nInp; j++) iDat[j] = inp(i,j);
            xc_eval(this->functional, iDat, oDat);
            for (int j = 0; j < nOut; j++) out(i,j) = oDat[j];
        } else {
            for (int j = 0; j < nOut; j++) out(i,j) = 0.0;
        }
    }
    delete[] iDat;
    delete[] oDat;
}


/*
  int xc_eval_setup(xc_functional fun,
		    enum xc_vars vars,
		    enum xc_mode mode,
		    int order);

  void xc_eval(xc_functional fun,
	       const double *density,
	       double *result);

  void xc_eval_vec(xc_functional fun, int nr_points,
		   const double *density,
		   int density_pitch,
		   double *result,
		   int result_pitch);

*/
