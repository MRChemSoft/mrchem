#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/trees/FunctionNode.h"

#include "XCFunctional.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionNode;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::DerivativeOperator;

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrdft {

/** @brief Constructor
 *
 * Initializes the new functional
 *
 * @param[in] mra Computational domain for the density grid
 * @param[in] spin True for spin-separated calculations
 * @param[in] order Order of XC kernel derivatives
 *
 */
XCFunctional::XCFunctional(mrcpp::MultiResolutionAnalysis<3> &mra, bool spin)
        : spin_separated(spin),
          MRA(mra),
          use_gamma(false),
          cutoff(-1.0),
          functional(xc_new_functional()),
          derivative(nullptr),
          rho_a(nullptr),
          rho_b(nullptr),
          rho_t(nullptr) {
    derivative = new mrcpp::ABGVOperator<3>(this->MRA, 0.0, 0.0);
    if (isSpinSeparated()) {
        rho_a = new FunctionTree<3>(this->MRA);
        rho_b = new FunctionTree<3>(this->MRA);
    } else {
        rho_t = new FunctionTree<3>(this->MRA);
    }
}

/** @brief Destructor */
XCFunctional::~XCFunctional() {
    xc_free_functional(functional);
    if (rho_a != nullptr) delete rho_a;
    if (rho_b != nullptr) delete rho_b;
    if (rho_t != nullptr) delete rho_t;
    if (derivative != nullptr) delete derivative;
}

/** @brief User-friendly setup of the xcfun calculation
 *
 * Setup the XC functional for evaluation. In MRChem we use only a subset
 * of the alternatives offered by xcfun. More functionality might be enabled
 * at a later stage.
 *
 * @param[in] order Order of the requested operator (1 for potential, 2 for hessian, ...)
 *
 */
void XCFunctional::evalSetup(int order) {
    unsigned int func_type = isGGA();               //!< only LDA and GGA supported for now
    unsigned int dens_type = 1 + isSpinSeparated(); //!< only n (dens_type = 1) or alpha & beta (denst_type = 2) supported now.
    unsigned int mode_type = 1;                     //!< only derivatives (neither potential nor contracted)
    unsigned int laplacian = 0;                     //!< no laplacian
    unsigned int kinetic = 0;                       //!< no kinetic energy density
    unsigned int current = 0;                       //!< no current density
    unsigned int exp_derivative = not useGamma();   //!< use gamma or explicit derivatives
    if (isLDA()) exp_derivative = 0;                //!< fall back to gamma-type derivatives if LDA (bad hack: no der are actually needed here!)
    xc_user_eval_setup(functional, order, func_type, dens_type, mode_type, laplacian, kinetic, current, exp_derivative);
}

/** @brief Functional setup
 *
 * @param[in] name The name of the chosen functional
 * @param[in] coef The amount of the chosen functional
 *
 * For each functional part in calculation a corresponding token is
 * created in xcfun.
 *
 */
void XCFunctional::setFunctional(const std::string &name, double coef) {
    xc_set(functional, name.c_str(), coef);
}

bool XCFunctional::hasDensity() const {
    bool out = true;
    if (isSpinSeparated()) {
        if (rho_a->getSquareNorm() < 0.0) out = false;
        if (rho_b->getSquareNorm() < 0.0) out = false;
    } else {
        if (rho_t->getSquareNorm() < 0.0) out = false;
    }
    return out;
}

/** @brief Return FunctionTree for the xcfun input density
 *
 * @param[in] type Which density to return (alpha, beta or total)
 *
 * Returns a reference to the internal density function so that
 * it can be projected by the host program. This needs to be done
 * before setup().
 */
FunctionTree<3> & XCFunctional::getDensity(DensityType type) {
    switch (type) {
    case DensityType::Total:
        if (rho_t == nullptr) MSG_FATAL("Total density not allocated");
        return *rho_t;
    case DensityType::Alpha:
        if (rho_a == nullptr) MSG_FATAL("Alpha density not allocated");
        return *rho_a;
    case DensityType::Beta:
        if (rho_b == nullptr) MSG_FATAL("Beta density not allocated");
        return *rho_b;
    default:
        MSG_FATAL("Invalid density type");
        break;
    }
}

/** @brief Prepare for xcfun evaluation
 *
 * This computes the necessary input functions (gradients and gamma) and
 * constructs empty grids to hold the output functions of xcfun, and
 * collects them in the xcInput and xcOutput vectors, respectively.
 * Assumes that the density functions rho_t or rho_a/rho_b are already computed.
 */
void XCFunctional::setup() {
    if (isGGA()) {
        if (isSpinSeparated()) {
            grad_a = calcGradient(*rho_a);
            grad_b = calcGradient(*rho_b);
        } else {
            grad_t = calcGradient(*rho_t);
        }
    }
    if (useGamma()) {
        if (isSpinSeparated()) {
            gamma.push_back(calcDotProduct(grad_a, grad_a));
            gamma.push_back(calcDotProduct(grad_a, grad_b));
            gamma.push_back(calcDotProduct(grad_b, grad_b));
        } else {
            gamma.push_back(calcDotProduct(grad_t, grad_t));
        }
    }
    setupXCInput();
    setupXCOutput();
}

/** @brief Cleanup after xcfun evaluation
 *
 * This deallocates the memory used by both the input and output functions
 * to xcfun. The output must be collected (calcEnergy() and calcPotential())
 * before this function is called.
 */
void XCFunctional::clear() {
    xcInput.clear(false);
    xcOutput.clear(true);
    grad_a.clear(true);
    grad_b.clear(true);
    grad_t.clear(true);
    gamma.clear(true);
    if (rho_a != nullptr) rho_a->clear(); // resets the tree to only root nodes
    if (rho_b != nullptr) rho_b->clear(); // resets the tree to only root nodes
    if (rho_t != nullptr) rho_t->clear(); // resets the tree to only root nodes
}

/** @brief Allocate input arrays for xcfun
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
    int nUsed = setupXCInputDensity();
    if (isGGA()) nUsed += setupXCInputGradient();

    if (nInp != nUsed) MSG_ERROR("Mismatch between used vs requested " << nUsed << " : " << nInp);
    for (int i = 0; i < nInp; i++) {
        if (xcInput[i] == 0) MSG_ERROR("Invalid XC input");
    }
}

/** @brief Sets xcInput pointers for the density
 *
 * Returns the nr. of pointers used for sanity checking
 */
int XCFunctional::setupXCInputDensity() {
    int nUsed = 0;
    if (isSpinSeparated()) {
        if (rho_a == nullptr) MSG_FATAL("Invalid alpha density");
        if (rho_b == nullptr) MSG_FATAL("Invalid beta density");
        xcInput.push_back(rho_a);
        xcInput.push_back(rho_b);
        nUsed = 2;
    } else {
        if (rho_t == nullptr) MSG_FATAL("Invalid total density");
        xcInput.push_back(rho_t);
        nUsed = 1;
    }
    return nUsed;
}

/** @brief Sets xcInput pointers for the gradient(s)
 *
 * Returns the nr. of pointers used for sanity checking
 */
int XCFunctional::setupXCInputGradient() {
    int nUsed = 0;
    if (useGamma()) {
        xcInput.push_back(gamma);
        nUsed += gamma.size();
    } else {
        if (isSpinSeparated()) {
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

/** @brief Allocate output arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocated and the function objects are created, borrowing the
 * grid from the electronic density.
 *
 */
void XCFunctional::setupXCOutput() {
    if (xcOutput.size() != 0) MSG_ERROR("XC output not empty");
    if (xcInput.size() == 0) MSG_ERROR("XC input not initialized");

    // Fetch input density to copy its grid
    FunctionTree<3> &rho = *xcInput[0];

    int nOut = getOutputLength();
    for (int i = 0; i < nOut; i++) {
        FunctionTree<3> * tmp = new FunctionTree<3>(this->MRA);
        mrcpp::copy_grid(*tmp, rho);
        xcOutput.push_back(tmp);
    }
}

/** @brief Evaluation of the functional and its derivatives
 *
 * The data contained in the xcInput is converted in matrix form and fed to
 * the functional. The output matrix is then converted back to function form.
 *
 */
void XCFunctional::evaluate() {
    if (xcInput.size() == 0) MSG_ERROR("XC input not initialized");
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");

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
            evaluateBlock(inpData, outData);
            expandNodeData(n, nOut, xcOutput, outData);
        }
    }
    for (int i = 0; i < nOut; i++) {
        xcOutput[i]->mwTransform(mrcpp::BottomUp);
        xcOutput[i]->calcSquareNorm();
    }

    timer.stop();
    Printer::printTree(0, "XC evaluate xcfun", xcOutput.sumNodes(), timer.getWallTime());
    printout(2, std::endl);
}

/** \brief Evaluates XC functional and derivatives
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
void XCFunctional::evaluateBlock(MatrixXd &inp, MatrixXd &out) const {
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

/** @brief Converts data from a FunctionNode to a matrix
 *
 * The FunctionNode(s) row data is packed into a matrix whose
 * dimensions are the overall number of grid points (nCoefs) and the
 * number of functions (nFuncs).
 *
 * parma[in] n the index of the requested node
 * param[in] nFuncs the number of functions
 * param[in] trees the array of FunctionTree(s)
 * param[in] the matrix object
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

/** @brief Converts data from a matrix to a FunctionNode
 *
 * The matrix containing the output from xcfun is converted back to the corresponding FunctionNode(s). The matrix dimensions are the overall number of grid points (nCoefs) and the number of functions (nFuncs).
 *
 * parma[in] n the index of the requested node
 * param[in] nFuncs the number of functions
 * param[in] trees the array of FunctionTree(s)
 * param[in] the matrix object
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

/** @brief Computes the XC energy as the integral of the functional */
double XCFunctional::calcEnergy() {
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");
    if (xcOutput[0] == 0) MSG_ERROR("Invalid XC output");

    Timer timer;
    double energy = xcOutput[0]->integrate();
    timer.stop();
    Printer::printTree(0, "XC energy", xcOutput[0]->getNNodes(), timer.getWallTime());
    return energy;
}

/** @brief Compute the XC potential
 *
 * Combines the xcfun output functions into the final XC potential functions.
 * Different calls for LDA and GGA.
 *
 */
FunctionTreeVector<3> XCFunctional::calcPotential() {
    FunctionTreeVector<3> xc_pot;
    if (xcOutput.size() == 0) MSG_ERROR("XC output not initialized");

    if (isLDA()) {
        calcPotentialLDA(xc_pot);
    } else if (isGGA()) {
        calcPotentialGGA(xc_pot);
    } else {
        MSG_FATAL("Invalid functional type");
    }
    return xc_pot;
}

/** @brief Potential calculation for LDA functionals
 *
 * The potential conicides with the xcfun output, which is then
 * assigned to the corresponding potential functions.
 *
 */
void XCFunctional::calcPotentialLDA(FunctionTreeVector<3> &potentials) {
    int order = 1; //HACK order needs to be passed!
    int nPotentials = isSpinSeparated() ? order + 1 : 1;

    if (order != 1) NOT_IMPLEMENTED_ABORT;

    for (int i = 0; i < nPotentials; i++) {
        int outputIndex = i + 1;
        if (xcOutput[outputIndex] == 0) MSG_ERROR("Invalid XC output");
        // Deep copy
        mrcpp::FunctionTreeVector<3> vec;
        vec.push_back(xcOutput[outputIndex]);
        FunctionTree<3> *pot = new FunctionTree<3>(this->MRA);
        mrcpp::copy_grid(*pot, vec);
        mrcpp::add(-1.0, *pot, vec);
        potentials.push_back(pot);
    }
}

/** @brief  potential calculation for GGA functionals
 *
 * the potential functions are assembled from the xcfun output functions
 * The method used depends on whether the functional is spin-separated
 * and whether explicit or gamma-type derivatives have been used in xcfun.
 *
 */
void XCFunctional::calcPotentialGGA(FunctionTreeVector<3> &potentials) {
    FunctionTree<3> * pot;
    if (isSpinSeparated()) {
        FunctionTree<3> & df_da = *xcOutput[1];
        FunctionTree<3> & df_db = *xcOutput[2];
        if (useGamma()) {
            FunctionTree<3> & df_dgaa = *xcOutput[3];
            FunctionTree<3> & df_dgab = *xcOutput[4];
            FunctionTree<3> & df_dgbb = *xcOutput[5];
            pot = calcPotentialGGA(df_da, df_dgaa, df_dgab, grad_a, grad_b);
            potentials.push_back(pot);
            pot = calcPotentialGGA(df_db, df_dgbb, df_dgab, grad_b, grad_a);
            potentials.push_back(pot);
        } else {
            FunctionTreeVector<3> df_dga;
            FunctionTreeVector<3> df_dgb;
            df_dga.push_back(xcOutput[3]);
            df_dga.push_back(xcOutput[4]);
            df_dga.push_back(xcOutput[5]);
            df_dgb.push_back(xcOutput[6]);
            df_dgb.push_back(xcOutput[7]);
            df_dgb.push_back(xcOutput[8]);
            pot = calcPotentialGGA(df_da, df_dga);
            potentials.push_back(pot);
            pot = calcPotentialGGA(df_db, df_dgb);
            potentials.push_back(pot);
        }
    } else {
        FunctionTree<3> & df_dt = *xcOutput[1];
        if (useGamma()) {
            FunctionTree<3> & df_dgamma = *xcOutput[2];
            pot = calcPotentialGGA(df_dt, df_dgamma, grad_t);
            potentials.push_back(pot);
        } else {
            FunctionTreeVector<3> df_dgt;
            df_dgt.push_back(xcOutput[2]);
            df_dgt.push_back(xcOutput[3]);
            df_dgt.push_back(xcOutput[4]);
            pot = calcPotentialGGA(df_dt, df_dgt);
            potentials.push_back(pot);
        }
    }
    pot = 0;
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
                                                 FunctionTreeVector<3> grad_rho) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(std::make_tuple(1.0, &df_drho));

    FunctionTree<3> *tmp = 0;
    tmp = calcGradDotPotDensVec(df_dgamma, grad_rho);
    funcs.push_back(std::make_tuple(-2.0, tmp));

    FunctionTree<3> * V = addPotentialContributions(funcs);
    mrcpp::clear(funcs, false);
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
                                                 FunctionTreeVector<3> grad_rhob) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(std::make_tuple(1.0, &df_drhoa));

    FunctionTree<3> *tmp1 = 0;
    tmp1 = calcGradDotPotDensVec(df_dgaa, grad_rhoa);
    funcs.push_back(std::make_tuple(-2.0, tmp1));

    FunctionTree<3> *tmp2 = 0;
    tmp2 = calcGradDotPotDensVec(df_dgab, grad_rhob);
    funcs.push_back(std::make_tuple(-1.0, tmp2));

    FunctionTree<3> * V = addPotentialContributions(funcs);
    mrcpp::clear(funcs, false);
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
                                                 FunctionTreeVector<3> & df_dgr) {
    FunctionTreeVector<3> funcs;
    funcs.push_back(std::make_tuple(1.0, &df_drho));

    FunctionTree<3> * tmp = calcDivergence(df_dgr);
    funcs.push_back(std::make_tuple(-1.0, tmp));

    FunctionTree<3> * V = addPotentialContributions(funcs);
    mrcpp::clear(funcs, false);
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
    mrcpp::build_grid(*V, contributions);
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
FunctionTree<3>* XCFunctional::calcDivergence(FunctionTreeVector<3> &inp) {
    if (derivative == 0) MSG_ERROR("No derivative operator");

    FunctionTreeVector<3> tmp_vec;
    for (int d = 0; d < 3; d++) {
        FunctionTree<3> *out_d = new FunctionTree<3>(MRA);
        mrcpp::apply(*out_d, *derivative, mrcpp::get_func(inp, d), d);
        tmp_vec.push_back(std::make_tuple(1.0, out_d));
    }
    FunctionTree<3> *out = new FunctionTree<3>(MRA);
    mrcpp::build_grid(*out, tmp_vec);
    mrcpp::add(-1.0, *out, tmp_vec, 0); // Addition on union grid
    mrcpp::clear(tmp_vec, true);
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
                                                     FunctionTreeVector<3> &rho) {
    FunctionTreeVector<3> vec;
    for (int d = 0; d < 3; d++) {
        Timer timer;
        FunctionTree<3> &rho_d = mrcpp::get_func(rho, d);
        FunctionTree<3> *Vrho = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*Vrho, rho_d);
        mrcpp::multiply(-1.0, *Vrho, 1.0, V, rho_d, -1);
        vec.push_back(std::make_tuple(1.0, Vrho));

        timer.stop();
        double t = timer.getWallTime();
        int n = Vrho->getNNodes();
        Printer::printTree(2, "Multiply", n, t);
    }

    Timer timer;
    FunctionTree<3> *result = calcDivergence(vec);
    mrcpp::clear(vec, true);

    timer.stop();
    double t = timer.getWallTime();
    int n = result->getNNodes();
    Printer::printTree(2, "Gradient", n, t);
    return result;
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
        FunctionTree<3> *gradient_comp = new FunctionTree<3>(this->MRA);
        mrcpp::apply(*gradient_comp, *derivative, function, d);
        gradient.push_back(std::make_tuple(1.0, gradient_comp));
    }
    return gradient;
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
        FunctionTree<3> &tree_a = mrcpp::get_func(vec_a, d);
        FunctionTree<3> &tree_b = mrcpp::get_func(vec_b, d);
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        mrcpp::build_grid(*out_d, tree_a);
        mrcpp::build_grid(*out_d, tree_b);
        mrcpp::multiply(-1.0, *out_d, 1.0, tree_a, tree_b, 0);
        out_vec.push_back(std::make_tuple(1.0, out_d));
    }
    FunctionTree<3> *out = new FunctionTree<3>(*MRA);
    build_grid(*out, out_vec);
    mrcpp::add(-1.0, *out, out_vec, -1);
    mrcpp::clear(out_vec, true);
    return out;
}

} //namespace mrdft
