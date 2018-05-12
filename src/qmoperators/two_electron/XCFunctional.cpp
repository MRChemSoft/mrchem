#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/trees/FunctionNode.h"

#include "XCFunctional.h"
#include "Orbital.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionNode;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::DerivativeOperator;

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA;

/** @brief constructor
 *
 * Initializes the new functional
 *
 * @param[in] s True for spin-separated calculations
 * @param[in] thrs Threshold for func calculation
 *
 */
XCFunctional::XCFunctional(bool s, bool e, double thrs, OrbitalVector &phi, DerivativeOperator<3> *D)
        : density(s, false), //HACK: shared set as false for now...
          derivative(D),
          spin(s),
          expDerivatives(0),
          cutoff(thrs),
          energy(0.0),
          functional(xc_new_functional()),
          orbitals(&phi) {
    if (e) expDerivatives = 1;
}

/** @brief destructor
 *
 */
XCFunctional::~XCFunctional() {
    xc_free_functional(functional);
}

/** @brief evaluation of the functional and its derivatives
 *
 * the data contained in the xcInput is converted in matrix form and fed to the functional.
 * the output matrix is then converted back to function form.
 *
 */
void XCFunctional::setup(const int order) {
    density::calc_density(density, *orbitals);
    if (this->isGGA()) calcDensityGradient(); // HACK we should implement gradient stuff as density-related functions, not here!
    if (this->needsGamma()) calcGamma();
    
    evalSetup(order);
    setupXCInput();
    setupXCOutput();
    potentialFunction.clear(true);
    evaluateFunctional();
    calcEnergy();
    calcPotential();
    density.clear();
    xcInput.clear(false);
    xcOutput.clear(false);
    grad_a.clear(true);
    grad_b.clear(true);
    grad_t.clear(true);
    gamma.clear(true);
}


/** @brief functional setup
 *
 * @usage For each functional part in calculation a corresponding token is created in xcfun
 *
 * @param[in] name The name of the chosen functional
 * @param[in] coef The amount of the chosen functional
 */
void XCFunctional::setFunctional(const std::string &name, double coef) {
    xc_set(functional, name.c_str(), coef);
}

/** @brief User-friendly setup of the xcfun calculation
 *
 * Setup the XC functional for evaluation. In MRChem we use only a subset of the alternatives offered by xcfun.
 * More functinality might be enabled at a later stage.
 *
 * @param[in] order Order of the requested operator (1 for potential, 2 for hessian, ...)
 *
 */
void XCFunctional::evalSetup(const int order) {
    unsigned int func_type = this->isGGA();  //!< only LDA and GGA supported for now
    unsigned int dens_type = 1 + this->isSpinSeparated(); //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int mode_type = 1; //!< only derivatives (neither potential nor contracted)
    unsigned int laplacian = 0; //!< no laplacian
    unsigned int kinetic = 0;   //!< no kinetic energy density
    unsigned int current = 0;   //!< no current density
    if(this->isLDA()) { // Fall back to gamma-type derivatives if LDA (bad hack: no der are actually needed here!)
        expDerivatives = 0;
    }
    xc_user_eval_setup(functional, order, func_type, dens_type, mode_type, laplacian, kinetic, current, expDerivatives);
}


/** \breif Evaluates XC functional and derivatives
 *
 * Computes the alpha and beta exchange-correlation functionals and
 * their derivatives.  The electronic density (total/alpha/beta) and their gradients are
 * given as input. Results are then stored in the xcfun output
 * functions. Higher order derivatives can be computed changing the parameter k.
 *
 * XCFunctional output (with k=1 and explicit derivatives):
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
 *
 * XCFunctional output (with k=1 and gamma-type derivatives):
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \gamma} \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta },
 *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \alpha}},
 *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \beta }},
 *  \frac{\partial F_{xc}}{\partial \gamma^{\beta  \beta }}
 *  \right) \f$
 *
 * The points are passed with through a matrix of dimension nInp x nPts
 * where nInp is the number of input data required for a single evaluation
 * and nPts is the number of points requested. Similarly the output is provided
 * as a matrix nOut x nPts.
 *
 * param[in] k the order of the requested derivatives
 * param[in] input values
 * param[out] output values
 *
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
        if (inp(i,0) > cutoff) {
            for (int j = 0; j < nInp; j++) iDat[j] = inp(i,j);
            xc_eval(functional, iDat, oDat);
            for (int j = 0; j < nOut; j++) out(i,j) = oDat[j];
        } else {
            for (int j = 0; j < nOut; j++) out(i,j) = 0.0;
        }
    }
    delete[] iDat;
    delete[] oDat;
}

/** @brief XC potential calculation
 *
 * Computes the XC potential for a non-spin separated functional and
 * gamma-type derivatives
 *
 * @param[in] df_drho functional derivative wrt rho
 * @param[in] df_dgamma functional_derivative wrt gamma
 * @param[in] grad_rho gradient of rho
 * @param[in] derivative derivative operator to use
 *
 */
FunctionTree<3> * XCFunctional::calcPotentialGGA(FunctionTree<3> & df_drho,
                                                 FunctionTree<3> & df_dgamma,
                                                 FunctionTreeVector<3> grad_rho,
                                                 DerivativeOperator<3> *derivative) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, &df_drho);

    FunctionTree<3> *tmp = 0;
    tmp = calcGradDotPotDensVec(df_dgamma, grad_rho, derivative);
    funcs.push_back(-2.0, tmp);

    FunctionTree<3> * V = addPotentialContributions(funcs);
    funcs.clear(false);
    delete tmp;
    return V;
}

/** @brief XC potential calculation
 *
 * Computes the XC potential for a spin separated functional and
 * gamma-type derivatives
 *
 * @param[in] df_drhoa functional derivative wrt rhoa
 * @param[in] df_drhob functional derivative wrt rhob
 * @param[in] df_dgaa  functional_derivative wrt gamma_aa
 * @param[in] df_dgab  functional_derivative wrt gamma_ab
 * @param[in] df_dgbb  functional_derivative wrt gamma_bb
 * @param[in] grad_rhoa gradient of rho_a
 * @param[in] grad_rhob gradient of rho_b
 * @param[in] derivative derivative operator to use
 *
 */
FunctionTree<3> * XCFunctional::calcPotentialGGA(FunctionTree<3> & df_drhoa,
                                                 FunctionTree<3> & df_dgaa,
                                                 FunctionTree<3> & df_dgab,
                                                 FunctionTreeVector<3> grad_rhoa,
                                                 FunctionTreeVector<3> grad_rhob,
                                                 DerivativeOperator<3> *derivative) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, &df_drhoa);

    FunctionTree<3> *tmp1 = 0;
    tmp1 = calcGradDotPotDensVec(df_dgaa, grad_rhoa, derivative);
    funcs.push_back(-2.0, tmp1);

    FunctionTree<3> *tmp2 = 0;
    tmp2 = calcGradDotPotDensVec(df_dgab, grad_rhob, derivative);
    funcs.push_back(-1.0, tmp2);

    FunctionTree<3> * V = addPotentialContributions(funcs);
    funcs.clear(false);
    delete tmp1;
    delete tmp2;
    return V;
}

/** @brief XC potential calculation
 *
 * Computes the XC potential for explicit derivatives.
 *
 * @param[in] df_drho functional derivative wrt rho
 * @param[in] df_dgr  functional_derivative wrt grad_rho
 * @param[in] derivative derivative operator to use
 *
 */
FunctionTree<3> * XCFunctional::calcPotentialGGA(FunctionTree<3> & df_drho,
                                                 FunctionTreeVector<3> & df_dgr,
                                                 DerivativeOperator<3> *derivative) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, &df_drho);

    FunctionTree<3> * tmp = calcDivergence(df_dgr, derivative);
    funcs.push_back(-1.0, tmp);

    FunctionTree<3> * V = addPotentialContributions(funcs);
    funcs.clear(false);
    delete tmp;
    return V;
}

/** @brief adds all potential contributions together
 *
 * @param[in] contributions vctor with all contributions
 *
 * NOTE: this should possibly be moved to the new mrcpp module
 * as it only involves mwtrees
 */
FunctionTree<3> * XCFunctional::addPotentialContributions(FunctionTreeVector<3> & contributions) {
    FunctionTree<3> *V = new FunctionTree<3>(*MRA);
    mrcpp::copy_grid(*V, contributions);
    mrcpp::add(-1, *V, contributions, -1);
    return V;
}

/** @brief computes the divergence of a vector field
 *
 * @param[in] inp the vector field expressed as function trees
 * @param[in] derivative the derivative operator
 *
 * NOTE: this should possibly be moved to the new mrcpp module
 * as it only involves mwtrees
 */
FunctionTree<3>* XCFunctional::calcDivergence(FunctionTreeVector<3> &inp,
                                              DerivativeOperator<3> *derivative) {
    if (derivative == 0) MSG_ERROR("No derivative operator");

    FunctionTreeVector<3> tmp_vec;
    for (int d = 0; d < 3; d++) {
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        mrcpp::apply(*out_d, *derivative, *inp[d], d);
        tmp_vec.push_back(out_d);
    }
    FunctionTree<3> *out = new FunctionTree<3>(*MRA);
    mrcpp::copy_grid(*out, tmp_vec);
    mrcpp::add(-1.0, *out, tmp_vec, 0); // Addition on union grid
    tmp_vec.clear(true);
    return out;
}

/** brief divergenge of a vector field times a function
 *
 * @param[in]V Function (derivative of the functional wrt gamma)
 * @param[in]rho vector field (density gradient)
 * @param[in] derivative the derivative operator
 *
 * NOTE: this should possibly be moved to the new mrcpp module
 * as it only involves mwtrees
 */
FunctionTree<3>* XCFunctional::calcGradDotPotDensVec(FunctionTree<3> &V,
                                                     FunctionTreeVector<3> &rho,
                                                     DerivativeOperator<3> *derivative) {
    FunctionTreeVector<3> vec;
    for (int d = 0; d < 3; d++) {
        if (rho[d] == 0) MSG_ERROR("Invalid density");

        Timer timer;
        FunctionTree<3> *Vrho = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*Vrho, *rho[d]);
        mrcpp::multiply(-1.0, *Vrho, 1.0, V, *rho[d], -1);
        vec.push_back(Vrho);

        timer.stop();
        double t = timer.getWallTime();
        int n = Vrho->getNNodes();
        Printer::printTree(2, "Multiply", n, t);
    }

    Timer timer;
    FunctionTree<3> *result = calcDivergence(vec, derivative);
    vec.clear(true);

    timer.stop();
    double t = timer.getWallTime();
    int n = result->getNNodes();
    Printer::printTree(2, "Gradient", n, t);
    return result;
}

/** @brief allocate input arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocared and its pointers assigned to the required input
 * functions.
 *
 */
void XCFunctional::setupXCInput() {
    if (xcInput.size() != 0) MSG_ERROR("XC input not empty");
    println(2, "Preprocessing");
    
    int nInp = getInputLength();
    
    int nUsed = 0;
    nUsed = setupXCInputDensity(nUsed);
    if (this->isGGA()) {
        nUsed = setupXCInputGradient(nUsed);
    }
    if (nInp != nUsed)  MSG_ERROR("Mismatch between used vs requested");
    for (int i = 0; i < nInp; i++) {
        if (xcInput[i] == 0) MSG_ERROR("Invalid XC input");
    }
}

/** @brief sets xcInput pointers for the density
 *
 * Returns the nr. of pointers used for sanity checking
 */
int XCFunctional::setupXCInputDensity(int nUsed) {
    if (this->isSpinSeparated()) {
        xcInput.push_back(&density.alpha());
        xcInput.push_back(&density.beta());
        nUsed += 2;
    } else {
        xcInput.push_back(&density.total());
        nUsed++;
    }
    return nUsed;
}

/** @brief sets xcInput pointers for the gradient(s)
 *
 * Returns the nr. of pointers used for sanity checking
 */
int XCFunctional::setupXCInputGradient(int nUsed) {
    if(this->needsGamma()) {
        xcInput.push_back(gamma);
        nUsed += gamma.size();
    } else {
        if(this->isSpinSeparated()) {
            xcInput.push_back(grad_a);
            xcInput.push_back(grad_b);
            nUsed += grad_a.size() + grad_b.size();
        } else {
            xcInput.push_back(grad_t);
            nUsed += grad_t.size();
        }
    }
    return nUsed;
}

/** @brief allocate output arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocated and the function objects are created, borrowing the
 * grid from the electronic density.
 *
 */
void XCFunctional::setupXCOutput() {
    if (xcOutput.size() != 0) MSG_ERROR("XC output not empty");
    if (xcInput.size() == 0) MSG_ERROR("XC input not initialized");

    // Alloc output trees
    int nOut = getOutputLength();

    // Copy grid from input density
    FunctionTree<3> &rho = *xcInput[0];
    for (int i = 0; i < nOut; i++) {
        FunctionTree<3> * tmp = new FunctionTree<3>(*MRA);
        copy_grid(*tmp, rho);
        xcOutput.push_back(tmp);
    }
}

/** @brief evaluation of the functional and its derivatives
 *
 * the data contained in the xcInput is converted in matrix form and fed to the functional.
 * the output matrix is then converted back to function form.
 *
 */
void XCFunctional::evaluateFunctional() {
    if (xcInput.size() == 0) MSG_ERROR("XC input not initialized");
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");

    int order = 1; //HACK order needs to be passed!

    Timer timer;
    println(2, "Evaluating");

    int nInp = getInputLength();
    int nOut = getOutputLength();

#pragma omp parallel firstprivate(nInp, nOut)
    {
        int nNodes = xcInput[0]->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MatrixXd inpData, outData;
            compressNodeData(n, nInp, xcInput, inpData);
            evaluate(order, inpData, outData);
            expandNodeData(n, nOut, xcOutput, outData);
        }
    }
    for (int i = 0; i < nOut; i++) {
        xcOutput[i]->mwTransform(mrcpp::BottomUp);
        xcOutput[i]->calcSquareNorm();
    }

    timer.stop();
    double t = timer.getWallTime();
    int n = xcOutput.sumNodes();
    Printer::printTree(0, "XC evaluate xcfun", n, t);
    printout(2, std::endl);
}

/** @brief converts data from a FunctionNode to a matrix
 *
 * The FunctionNode(s) row data is packed into a matrix whose
 * dimensions are the overall number of grid points (nCoefs) and the
 * number of functions (nFuncs).
 *
 * parma[in] n the index of the requested node
 * param[in] nFuncs the number of functions
 * param[in] trees the array of FunctionTree(s)
 * param[in] the matrix object.
 */
void XCFunctional::compressNodeData(int n, int nFuncs, FunctionTreeVector<3> trees, MatrixXd &data) {
    if (trees.size() == 0) MSG_ERROR("Invalid input");

    FunctionTree<3> &tree = *trees[0];
    int nCoefs = tree.getTDim()*tree.getKp1_d();
    data = MatrixXd::Zero(nCoefs, nFuncs);

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized input tree");
        FunctionNode<3> &node = trees[i]->getEndFuncNode(n);
        VectorXd col_i;
        node.getValues(col_i);
        data.col(i) = col_i;
    }
}

/** @brief converts data from a matrix to a FunctionNode
 *
 * The matrix containing the output from xcfun is converted back to the corresponding FunctionNode(s). The matrix dimensions are the overall number of grid points (nCoefs) and the number of functions (nFuncs).
 *
 * parma[in] n the index of the requested node
 * param[in] nFuncs the number of functions
 * param[in] trees the array of FunctionTree(s)
 * param[in] the matrix object.
 */
void XCFunctional::expandNodeData(int n, int nFuncs, FunctionTreeVector<3> trees, MatrixXd &data) {
    if (trees.size() == 0) MSG_ERROR("Invalid input");

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized output tree " << i);
        VectorXd col_i = data.col(i);
        FunctionNode<3> &node = trees[i]->getEndFuncNode(n);
        node.setValues(col_i);
    }
}

/** @brief potential calculation
 *
 * different calls for LDA and GGA
 *
 */
void XCFunctional::calcPotential() {
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");
    
    if (this->isLDA()) {
        calcPotentialLDA();
    } else if (this->isGGA()) {
        calcPotentialGGA();
    } else {
        MSG_FATAL("Invalid functional type");
    }
}

/** @brief potential calculation for LDA functionals
 *
 * The potential conicides with the xcfun output, which is then
 * assigned to the corresponding potential functions.
 *
 */
void XCFunctional::calcPotentialLDA() {
    int order = 1; //HACK order needs to be passed!
    int nPotentials = this->isSpinSeparated() ? order + 1 : 1;

    if (order != 1) NOT_IMPLEMENTED_ABORT;

    for (int i = 0; i < nPotentials; i++) {
        int outputIndex = i + 1;
        if (xcOutput[outputIndex] == 0) MSG_ERROR("Invalid XC output");
        potentialFunction.push_back(xcOutput[outputIndex]);
    }
}

/** @brief  potential calculation for GGA functionals
 *
 * the potential functions are assembled from the xcfun output functions
 * The method used depends on whether the functional is spin-separated
 * and whether explicit or gamma-type derivatives have been used in xcfun.
 *
 */
void XCFunctional::calcPotentialGGA() {
    FunctionTree<3> * pot;
    if (this->isSpinSeparated()) {
        FunctionTree<3> & df_da = *xcOutput[1];
        FunctionTree<3> & df_db = *xcOutput[2];
        if(this->needsGamma()) {
            FunctionTree<3> & df_dgaa = *xcOutput[3];
            FunctionTree<3> & df_dgab = *xcOutput[4];
            FunctionTree<3> & df_dgbb = *xcOutput[5];
            pot = calcPotentialGGA(df_da, df_dgaa, df_dgab, grad_a, grad_b, derivative);
            potentialFunction.push_back(pot);
            pot = calcPotentialGGA(df_db, df_dgbb, df_dgab, grad_b, grad_a, derivative);
            potentialFunction.push_back(pot);
        }
        else {
            FunctionTreeVector<3> df_dga;
            FunctionTreeVector<3> df_dgb;
            df_dga.push_back(xcOutput[3]);
            df_dga.push_back(xcOutput[4]);
            df_dga.push_back(xcOutput[5]);
            df_dgb.push_back(xcOutput[6]);
            df_dgb.push_back(xcOutput[7]);
            df_dgb.push_back(xcOutput[8]);
            pot = calcPotentialGGA(df_da, df_dga, derivative);
            potentialFunction.push_back(pot);
            pot = calcPotentialGGA(df_db, df_dgb, derivative);
            potentialFunction.push_back(pot);
        }

    }
    else {
        FunctionTree<3> & df_dt = *xcOutput[1];
        if(this->needsGamma()) {
            FunctionTree<3> & df_dgamma = *xcOutput[2];
            pot = calcPotentialGGA(df_dt, df_dgamma, grad_t, derivative);
            potentialFunction.push_back(pot);
        }
        else {
            FunctionTreeVector<3> df_dgt;
            df_dgt.push_back(xcOutput[2]);
            df_dgt.push_back(xcOutput[3]);
            df_dgt.push_back(xcOutput[4]);
            pot = calcPotentialGGA(df_dt, df_dgt, derivative);
            potentialFunction.push_back(pot);
        }
    }
    pot = 0;
}

/** @brief computes the gradient invariants (gamma functions)
 *
 * Depending on the mode chosen, xcfun needs either the gamma
 * functions or the explicit gradients. The first mode is possibly
 * more efficient (fewer functions to compute/handle), whereas the
 * other is simpler to implement. We keep both options open and
 * compute the gradient invariants if and when necessary.
 *
 */
void XCFunctional::calcGamma() {
    FunctionTree<3> * temp;
    if (this->isSpinSeparated()) {
        temp = calcDotProduct(grad_a, grad_a);
        gamma.push_back(temp);
        temp = calcDotProduct(grad_a, grad_b);
        gamma.push_back(temp);
        temp = calcDotProduct(grad_b, grad_b);
        gamma.push_back(temp);
    }
    else {
        temp = calcDotProduct(grad_t, grad_t);
        gamma.push_back(temp);
    }
}

/** @brief computes the gradient of the density
 *
 * For spin-free calculations, the total density is used. For
 * spin-separated calculations both alpha and beta gradients are
 * computed. The results are stored in the correspondig data members
 * of the XCFunctional.
 *
 */
int XCFunctional::calcDensityGradient() {
    int nNodes = 0;
    if (density.isSpinDensity()) {
        grad_a = calcGradient(density.alpha());
        grad_b = calcGradient(density.beta());
        nNodes += grad_a[0]->getNNodes() + grad_a[1]->getNNodes() + grad_a[2]->getNNodes();
        nNodes += grad_b[0]->getNNodes() + grad_b[1]->getNNodes() + grad_b[2]->getNNodes();
    }
    else {
        grad_t = calcGradient(density.total());
        nNodes += grad_t[0]->getNNodes() + grad_t[1]->getNNodes() + grad_t[2]->getNNodes();
    }
    return nNodes;
}

/** @brief computes and stores the gradient of a function
 *
 * @param[in] function
 *
 * Note: this should also be handled at a lower level (mrcpp)
 *
 */
FunctionTreeVector<3> XCFunctional::calcGradient(FunctionTree<3> &function) {
    if (derivative == 0) MSG_ERROR("No derivative operator");
    FunctionTreeVector<3> gradient;
    for (int d = 0; d < 3; d++) {
        FunctionTree<3> *gradient_comp = new FunctionTree<3>(*MRA);
        mrcpp::apply(*gradient_comp, *derivative, function, d);
        gradient.push_back(gradient_comp);
    }
    return gradient;
}

/** @brief computes the XC energy as the integral of the functional.
 *
 */
void XCFunctional::calcEnergy() {
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");
    if (xcOutput[0] == 0) MSG_ERROR("Invalid XC output");
    
    Timer timer;
    energy = xcOutput[0]->integrate();
    timer.stop();
    double t = timer.getWallTime();
    int n = xcOutput[0]->getNNodes();
    Printer::printTree(0, "XC energy", n, t);
}

/** @brief scalar product of two FunctionTreeVector(s)
 *
 * param[in] vec_a first vector
 * param[in] vec_b second vector
 *
 * Note: should be a mrcpp functionality.
 *
 */
FunctionTree<3>* XCFunctional::calcDotProduct(FunctionTreeVector<3> &vec_a,
                                              FunctionTreeVector<3> &vec_b) {
    if (vec_a.size() != vec_b.size()) MSG_ERROR("Invalid input");

    FunctionTreeVector<3> out_vec;
    for (int d = 0; d < vec_a.size(); d++) {
        FunctionTree<3> &tree_a = vec_a.getFunc(d);
        FunctionTree<3> &tree_b = vec_b.getFunc(d);
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*out_d, tree_a);
        mrcpp::copy_grid(*out_d, tree_b);
        mrcpp::multiply(-1.0, *out_d, 1.0, tree_a, tree_b, 0);
        out_vec.push_back(out_d);
    }
    FunctionTree<3> *out = new FunctionTree<3>(*MRA);
    copy_grid(*out, out_vec);
    mrcpp::add(-1.0, *out, out_vec, -1);

    out_vec.clear(true);
    return out;
}

/** @brief fetches the correct index for the potential function to use
 *
 * @param[in] orb the potentialFunction will be applied to this orbital.
 *
 * Based on the orbital spin, and whether the functional is spin
 * separated, the correct potential index is selected.
 *
 */
int XCFunctional::getPotentialFunctionIndex(const Orbital &orb) {
    int orbitalSpin = orb.spin();
    bool spinSeparatedFunctional = this->isSpinSeparated();
    int potentialFunctionIndex = -1;
    if (spinSeparatedFunctional and orbitalSpin == SPIN::Alpha) {
        potentialFunctionIndex = 0;
    }
    else if (spinSeparatedFunctional and orbitalSpin == SPIN::Beta) {
        potentialFunctionIndex = 1;
    }
    else if (not spinSeparatedFunctional) {
        potentialFunctionIndex = 0;
    }
    else {
        NOT_IMPLEMENTED_ABORT;
    }
    return potentialFunctionIndex;
}

} //namespace mrchem
