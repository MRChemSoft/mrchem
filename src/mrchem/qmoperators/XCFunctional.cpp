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
        xc_set_mode(this->functional, XC_VARS_AB);
    } else {
        xc_set_mode(this->functional, XC_VARS_N);
    }
}

XCFunctional::~XCFunctional() {
    xc_free_functional(this->functional);
}

void XCFunctional::setFunctional(const string &name, double coef) {
    int param = getParamFromName(name);
    xc_set_param(this->functional, param, coef);
}

/** \brief Obtains the parametrization of the functional based from
 * its name/keyword
 *
 */
int XCFunctional::getParamFromName(const string &name) {
    int param = -1;
    if (name == "SLATERX") {
        param = XC_SLATERX;
    } else if (name == "VWN5C") {
        param = XC_VWN5C;
    } else if (name == "BECKEX") {
        MSG_WARN("Functional not tested");
        param = XC_BECKEX;
    } else if (name == "BECKECORRX") {
        MSG_WARN("Functional not tested");
        param = XC_BECKECORRX;
    } else if (name == "BECKESRX") {
        MSG_WARN("Functional not tested");
        param = XC_BECKESRX;
    } else if (name == "OPTX") {
        MSG_WARN("Functional not tested");
        param = XC_OPTX;
    } else if (name == "LYPC") {
        MSG_WARN("Functional not tested");
        param = XC_LYPC;
    } else if (name == "PBEX") {
        param = XC_PBEX;
    } else if (name == "REVPBEX") {
        MSG_WARN("Functional not tested");
        param = XC_REVPBEX;
    } else if (name == "RPBEX") {
        MSG_WARN("Functional not tested");
        param = XC_RPBEX;
    } else if (name == "PBEC") {
        param = XC_PBEC;
    } else if (name == "SPBEC") {
        MSG_WARN("Functional not tested");
        param = XC_SPBEC;
    } else if (name == "VWN_PBEC") {
        MSG_WARN("Functional not tested");
        param = XC_VWN_PBEC;
    } else if (name == "LDAERFX") {
        MSG_WARN("Functional not tested");
        param = XC_LDAERFX;
    } else if (name == "LDAERFC") {
        MSG_WARN("Functional not tested");
        param = XC_LDAERFC;
    } else if (name == "LDAERFC_JT") {
        MSG_WARN("Functional not tested");
        param = XC_LDAERFC_JT;
    } else if (name == "KTX") {
        MSG_WARN("Functional not tested");
        param = XC_KTX;
    } else if (name == "TFK") {
        MSG_WARN("Functional not tested");
        param = XC_TFK;
    } else if (name == "PW91X") {
        MSG_WARN("Functional not tested");
        param = XC_PW91X;
    } else if (name == "PW91K") {
        MSG_WARN("Functional not tested");
        param = XC_PW91K;
    } else if (name == "PW92C") {
        MSG_WARN("Functional not tested");
        param = XC_PW92C;
    } else if (name == "MO5X") {
        MSG_WARN("Functional not tested");
        param = XC_M05X;
    } else if (name == "MO5X2X") {
        MSG_WARN("Functional not tested");
        param = XC_M05X2X;
    } else if (name == "MO6X") {
        MSG_WARN("Functional not tested");
        param = XC_M06X;
    } else if (name == "MO6X2X") {
        MSG_WARN("Functional not tested");
        param = XC_M06X2X;
    } else if (name == "MO6LX") {
        MSG_WARN("Functional not tested");
        param = XC_M06LX;
    } else if (name == "MO6HFX") {
        MSG_WARN("Functional not tested");
        param = XC_M06HFX;
    } else if (name == "BRX") {
        MSG_WARN("Functional not tested");
        param = XC_BRX;
    } else if (name == "MO5X2C") {
        MSG_WARN("Functional not tested");
        param = XC_M05X2C;
    } else if (name == "MO5C") {
        MSG_WARN("Functional not tested");
        param = XC_M05C;
    } else if (name == "MO6C") {
        MSG_WARN("Functional not tested");
        param = XC_M06C;
    } else if (name == "MO6HFC") {
        MSG_WARN("Functional not tested");
        param = XC_M06HFC;
    } else if (name == "MO6LC") {
        MSG_WARN("Functional not tested");
        param = XC_M06LC;
    } else if (name == "MO6X2C") {
        MSG_WARN("Functional not tested");
        param = XC_M06X2C;
    } else if (name == "TPSSC") {
        MSG_WARN("Functional not tested");
        param = XC_TPSSC;
    } else if (name == "TPSSX") {
        MSG_WARN("Functional not tested");
        param = XC_TPSSX;
    } else if (name == "REVTPSSC") {
        MSG_WARN("Functional not tested");
        param = XC_REVTPSSC;
    } else if (name == "REVTPSSX") {
        MSG_WARN("Functional not tested");
        param = XC_REVTPSSX;
    } else {
        MSG_ERROR("Invalid functional");
    }
    return param;
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
            xc_eval(this->functional, k, iDat, oDat);
            for (int j = 0; j < nOut; j++) out(i,j) = oDat[j];
        } else {
            for (int j = 0; j < nOut; j++) out(i,j) = 0.0;
        }
    }
    delete[] iDat;
    delete[] oDat;
}

