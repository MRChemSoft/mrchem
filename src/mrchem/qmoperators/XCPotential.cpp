#include "XCPotential.h"
#include "XCFunctional.h"
#include "QMPotential.h"
#include "FunctionTreeVector.h"
#include "GridGenerator.h"
#include "MWAdder.h"
#include "TelePrompter.h"
#include "Timer.h"

#include "constants.h"

extern MultiResolutionAnalysis<3> *MRA;

using namespace std;
using namespace Eigen;

/** \brief Creator
 *
 *  Initializes derivative order and index
 *
 */
XCPotential::XCPotential(int dI, int dO)
    : derivativeOrder(dO), derivativeIndex(dI) {
    if (abs(derivativeOrder) > 1) {
            NOT_IMPLEMENTED_ABORT;
        }
    if (abs(derivativeIndex) > abs(derivativeOrder)) {
        MSG_ERROR("Inconsistent input in XCPotential creator");
    }
	potentialFunction = 0;
}


void XCPotential::clear() {
    if(this->potentialFunction != 0) delete this->potentialFunction;
	this->potentialFunction = 0;
}


/** \brief Calculation of the GGA part of the XC potential
 *
 * The calculation is performed according the the following equation 
 * TO BE COMPLETED
 *
 * This function is identical irrespective of which case we consider
 * (alpha/beta/total). It is the caller which provides the right
 * assignment to the functions used here.
 */
FunctionTree<3>* XCPotential::calcPotentialGGA(FunctionTreeVector<3> &xc_funcs,
                                               FunctionTreeVector<3> &dRho_a,
                                               FunctionTreeVector<3> &dRho_b,
                                               DerivativeOperator<3> *derivative,
                                               int maxScale) {
    if (xc_funcs[0] == 0) MSG_ERROR("Invalid XC output");

    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, xc_funcs[0]);

    FunctionTree<3> *tmp_1 = 0;
    if (xc_funcs[1] != 0) {
        tmp_1 = calcGradDotPotDensVec(*xc_funcs[1], dRho_a, derivative, maxScale);
        funcs.push_back(-2.0, tmp_1);
    }

    FunctionTree<3> *tmp_2 = 0;
    if (xc_funcs[2] != 0) {
        tmp_2 = calcGradDotPotDensVec(*xc_funcs[2], dRho_b, derivative, maxScale);
        funcs.push_back(-1.0, tmp_2);
    }

    GridGenerator<3> G(maxScale);
    MWAdder<3> add(-1.0, maxScale);

    FunctionTree<3> *V = new FunctionTree<3>(*MRA);
    G(*V, funcs);
    add(*V, funcs, 0);
    funcs.clear(false);

    if (tmp_1 != 0) delete tmp_1;
    if (tmp_2 != 0) delete tmp_2;

    return V;
}

FunctionTree<3>* XCPotential::calcGradDotPotDensVec(FunctionTree<3> &V,
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

FunctionTree<3>* XCPotential::calcDivergence(FunctionTreeVector<3> &inp,
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

