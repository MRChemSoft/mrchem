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
    /*
    if (this->spin) {
        xc_set(this->functional, "LDA", 1.0);
    } else {
        xc_set(this->functional, "LDA", 1.0);
    }
    //LUCA: BAD HACK!!!
    xc_eval_setup(this->functional,
                  XC_N, XC_POTENTIAL,1);
    */
}

XCFunctional::~XCFunctional() {
    xc_free_functional(this->functional);
}

void XCFunctional::setFunctional(const string &name, double coef) {
    string funcName = "XC_" + name;
    xc_set(this->functional, funcName.c_str(), coef);
}

/*
  Setup the XC functional for evaluation. In MRChem we use only a subset of the alternatives offered by xcfun.
  More functinality might be enabled at a later stage.
 */
void XCFunctional::evalSetup(const int order)
{
    unsigned int func_type = this->isGGA(); //only LDA and GGA supported for now
    unsigned int dens_type = 1 + this->spin; // only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int mode_type = 2; // HACK POTENTIAL MODEcontracted mode only
    unsigned int laplacian = 0; // no laplacian
    unsigned int kinetic = 0; // no kinetic energy density
    unsigned int current = 0; // no current density
    unsigned int explicit_derivatives = 0; // only gamma-type derivatives now (soon to be changed!)
    std::cout << "In evalSetup. Order  " << order << " Dens " << dens_type << " Func " << func_type << std::endl;  
    xc_user_eval_setup(this->functional, order, func_type, dens_type, mode_type, laplacian, kinetic, current, explicit_derivatives);
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
    int nOut = getOutputLength(); // 2^order * n_points
    
    int nPts = inp.rows();
    out = MatrixXd::Zero(nPts, nOut);

    double *iDat = new double[nInp];
    double *oDat = new double[nOut];

    for (int i = 0; i < nPts; i++) {
        if (inp(i,0) > this->cutoff) {
            for (int j = 0; j < nInp; j++) iDat[j] = inp(i,j);
            xc_eval(this->functional, iDat, oDat);
            //std::cout << iDat[0] << oDat[0] << std::endl;
            for (int j = 0; j < nOut; j++) out(i,j) = oDat[j];
        } else {
            for (int j = 0; j < nOut; j++) out(i,j) = 0.0;
        }
    }
    delete[] iDat;
    delete[] oDat;
}
