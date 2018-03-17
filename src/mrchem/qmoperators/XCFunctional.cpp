#include "GridGenerator.h"
#include "MWAdder.h"
#include "MWMultiplier.h"
#include "MWDerivative.h"
#include "XCFunctional.h"
#include "TelePrompter.h"
#include "constants.h"

using namespace std;
using namespace Eigen;

extern MultiResolutionAnalysis<3> *MRA;

/** \brief creator 
 *
 * Initializes some parameters in the xcfun library
 *
 */
XCFunctional::XCFunctional(bool s, double thrs)
        : spin(s), cutoff(thrs) {
    this->functional = xc_new_functional();
    this->expDerivatives = 1; // explicit-type derivatives by default
}

XCFunctional::~XCFunctional() {
    xc_free_functional(this->functional);
}

void XCFunctional::setFunctional(const string &name, double coef) {
    xc_set(this->functional, name.c_str(), coef);
}

/*
  Setup the XC functional for evaluation. In MRChem we use only a subset of the alternatives offered by xcfun.
  More functinality might be enabled at a later stage.
 */
void XCFunctional::evalSetup(const int order)
{
    unsigned int func_type = this->isGGA(); //only LDA and GGA supported for now
    unsigned int dens_type = 1 + this->spin; // only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int mode_type = 1; //HACK HARD CODED
    unsigned int laplacian = 0; // no laplacian
    unsigned int kinetic = 0; // no kinetic energy density
    unsigned int current = 0; // no current density
    if(this->isLDA()) { // Fall back to gamma-type derivatives if LDA (bad hack: no der are actually needed here!)
        this->expDerivatives = 0;
    }
    xc_user_eval_setup(this->functional, order, func_type, dens_type, mode_type, laplacian, kinetic, current, this->expDerivatives);
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

FunctionTree<3> * XCFunctional::calcPotentialGGA(FunctionTree<3> & df_drho,
                                                 FunctionTree<3> & df_dgamma,
                                                 FunctionTreeVector<3> grad_rho,
                                                 DerivativeOperator<3> *derivative,
                                                 int maxScale) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, &df_drho);

    FunctionTree<3> *tmp = 0;
    tmp = calcGradDotPotDensVec(df_dgamma, grad_rho, derivative, maxScale);
    funcs.push_back(-2.0, tmp);

    FunctionTree<3> * V = addPotentialContributions(funcs, maxScale);
    funcs.clear(false);
    delete tmp;
    return V;
}

FunctionTree<3> * XCFunctional::calcPotentialGGA(FunctionTree<3> & df_drhoa,
                                                 FunctionTree<3> & df_dgaa,
                                                 FunctionTree<3> & df_dgab,
                                                 FunctionTreeVector<3> grad_rhoa,
                                                 FunctionTreeVector<3> grad_rhob,
                                                 DerivativeOperator<3> *derivative,
                                                 int maxScale) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, &df_drhoa);

    FunctionTree<3> *tmp1 = 0;
    tmp1 = calcGradDotPotDensVec(df_dgaa, grad_rhoa, derivative, maxScale);
    funcs.push_back(-2.0, tmp1);

    FunctionTree<3> *tmp2 = 0;
    tmp2 = calcGradDotPotDensVec(df_dgab, grad_rhob, derivative, maxScale);
    funcs.push_back(-1.0, tmp2);

    FunctionTree<3> * V = addPotentialContributions(funcs, maxScale);
    funcs.clear(false);
    delete tmp1;
    delete tmp2;
    return V;
}

FunctionTree<3> * XCFunctional::calcPotentialGGA(FunctionTree<3> & df_drho,
                                                 FunctionTreeVector<3> & df_dgr,
                                                 DerivativeOperator<3> *derivative,
                                                 int maxScale) {

    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, &df_drho);

    FunctionTree<3> * tmp = calcDivergence(df_dgr, derivative, maxScale);
    funcs.push_back(-1.0, tmp);

    FunctionTree<3> * V = addPotentialContributions(funcs, maxScale);
    funcs.clear(false);
    delete tmp;
    return V;
}

FunctionTree<3> * XCFunctional::addPotentialContributions(FunctionTreeVector<3> & contributions,
                                                         int maxScale) {
    FunctionTree<3> *V = new FunctionTree<3>(*MRA);
    GridGenerator<3> G(maxScale);
    MWAdder<3> add(-1.0, maxScale);
    G(*V, contributions);
    add(*V, contributions, 0);
    return V;
}

FunctionTree<3>* XCFunctional::calcDivergence(FunctionTreeVector<3> &inp,
                                             DerivativeOperator<3> *derivative,
                                             int maxScale) {
    if (derivative == 0) MSG_ERROR("No derivative operator");
    MWAdder<3> add(-1.0,  maxScale);
    MWDerivative<3> apply(maxScale);
    GridGenerator<3> grid(maxScale);

    FunctionTreeVector<3> tmp_vec;
    for (int d = 0; d < 3; d++) {
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        apply(*out_d, *derivative, *inp[d], d);
        tmp_vec.push_back(out_d);
    }
    FunctionTree<3> *out = new FunctionTree<3>(*MRA);
    grid(*out, tmp_vec);
    add(*out, tmp_vec, 0); // Addition on union grid
    tmp_vec.clear(true);
    return out;
}

FunctionTree<3>* XCFunctional::calcGradDotPotDensVec(FunctionTree<3> &V,
                                                    FunctionTreeVector<3> &rho,
                                                    DerivativeOperator<3> *derivative,
                                                    int maxScale) {
    MWMultiplier<3> mult(-1.0, maxScale);
    GridGenerator<3> grid(maxScale);

    FunctionTreeVector<3> vec;
    for (int d = 0; d < 3; d++) {
        if (rho[d] == 0) MSG_ERROR("Invalid density");

        Timer timer;
        FunctionTree<3> *Vrho = new FunctionTree<3>(*MRA);
        grid(*Vrho, *rho[d]);
        mult(*Vrho, 1.0, V, *rho[d], 0);
        vec.push_back(Vrho);

        timer.stop();
        double t = timer.getWallTime();
        int n = Vrho->getNNodes();
        TelePrompter::printTree(2, "Multiply", n, t);
    }

    Timer timer;
    FunctionTree<3> *result = calcDivergence(vec, derivative, maxScale);
    vec.clear(true);

    timer.stop();
    double t = timer.getWallTime();
    int n = result->getNNodes();
    TelePrompter::printTree(2, "Gradient", n, t);
    return result;
}

