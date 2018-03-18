#include "XCOperator.h"
#include "XCFunctional.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "Orbital.h"
#include "DensityProjector.h"
#include "Density.h"
#include "QMPotential.h"
#include "MWDerivative.h"
#include "TelePrompter.h"

using namespace std;
using namespace Eigen;

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param[in] k order of the operator
 * @param[in] F XCFunctional pointer
 * @param[in] phi vector of orbitals
 * @param[in] D derivative operators
 *
 * Based on the order and spin the correct nr. of potential functions is determined
 * Then the functional is set up for subsequent calculations, fixing some internals of
 * xcfun when F.evalSetup is invoked.
 *
 */
XCOperator::XCOperator(int k, XCFunctional &F, OrbitalVector &phi, DerivativeOperator<3> *D)
        : QMPotential(0),
          order(k),
          functional(&F),
          derivative(D),
          orbitals(&phi),
          energy(0.0),
          xcInput(0),
          xcOutput(0) {
    bool spin = F.isSpinSeparated();
    nPotentials = spin ? k + 1 : 1; /// k+1 potentials if spin separated, otherwise just one.
    density.setIsSpinDensity(spin);
    F.evalSetup(k);
}

/** @brief destructor
 *
 */
XCOperator::~XCOperator() {
    if (this->xcInput != 0)  MSG_ERROR("XC input not deallocated");
    if (this->xcOutput != 0) MSG_ERROR("XC output not deallocated");
    this->functional = 0;
    this->derivative = 0;
    this->orbitals = 0;
}

/** @brief setup the XCOperator
 * 
 * @param[in] prec precision 
 *
 * Sequence of steps required to compute the XC potentials The moat
 * important steps are evaluateXCFunctional and calcPotential where
 * the functional derivatives are computed and the potential assembled
 * respectively
 *
 */
void XCOperator::setup(double prec) {
    setApplyPrec(prec);
    calcDensity();
    setupXCInput();
    setupXCOutput();
    evaluateXCFunctional();
    calcEnergy();
    calcPotential();
    clearXCInput();
    clearXCOutput();
}

/** @brief XCOperator application
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
 *
 * @param[in] phi orbital to which the potential is applied.
 *
 */
Orbital* XCOperator::operator() (Orbital &phi) {
    FunctionTree<3> * potential = this->potentialFunction[this->getPotentialFunctionIndex(phi)];
    this->setReal(potential);
    this->setImag(0);
    Orbital * Vphi = QMPotential::operator()(phi); 
    this->clearReal(false);
    this->clearImag(false);
    return Vphi;
}

/** @brief adjoint operator not implemented
 *
 */
Orbital* XCOperator::adjoint(Orbital &phi) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief potential calculation
 *
 * different calls for LDA and GGA
 *
 */
void XCOperator::calcPotential() {
    if (xcOutput == 0) MSG_ERROR("XC output not initialized");
    
    bool lda = this->functional->isLDA();
    bool gga = this->functional->isGGA();
    
    if (lda) {
        calcPotentialLDA();
    } else if (gga) {
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
void XCOperator::calcPotentialLDA() {

    if (this->order != 1) {
        NOT_IMPLEMENTED_ABORT;
    }
    for (int i = 0; i < this->nPotentials; i++) {
        int outputIndex = i + 1;
        if (xcOutput[outputIndex] == 0) MSG_ERROR("Invalid XC output");
        potentialFunction.push_back(xcOutput[outputIndex]);
        xcOutput[outputIndex] = 0;
    }
}

/** @brief  potential calculation for GGA functionals
 *
 * the potential functions are assembled from the xcfun output functions
 * The metod used depends on whether the functional is spin-separated 
 * and whether explicit or gamma-type derivatives have been used in xcfun.
 * The corresponding method in the XCFunctional class is then selected and used
 *
 * Note: maybe all this stuff should end up in XCFunctional so that we
 * don't need to expose XCOutput.
 */
void XCOperator::calcPotentialGGA() {

    FunctionTree<3> * pot;
    bool spin = this->functional->isSpinSeparated();
    bool needsGamma = this->functional->needsGamma();
    if(spin) {
        FunctionTree<3> & df_da   = *xcOutput[1];
        FunctionTree<3> & df_db   = *xcOutput[2];
        if(needsGamma) {
            FunctionTree<3> & df_dgaa = *xcOutput[3];
            FunctionTree<3> & df_dgab = *xcOutput[4];
            FunctionTree<3> & df_dgbb = *xcOutput[5];
            pot = this->functional->calcPotentialGGA(df_da, df_dgaa, df_dgab, grad_a, grad_b,
                                                     this->derivative, this->max_scale);
            potentialFunction.push_back(pot);
            pot = this->functional->calcPotentialGGA(df_db, df_dgbb, df_dgab, grad_b, grad_a,
                                                     this->derivative, this->max_scale);
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
            pot = this->functional->calcPotentialGGA(df_da, df_dga, this->derivative,
                                                     this->max_scale);
            potentialFunction.push_back(pot);
            pot = this->functional->calcPotentialGGA(df_db, df_dgb, this->derivative,
                                                     this->max_scale);
            potentialFunction.push_back(pot);
        }
            
    }
    else {
        FunctionTree<3> & df_dt = *xcOutput[1];
        if(needsGamma) {
            FunctionTree<3> & df_dgamma = *xcOutput[2];
            pot = this->functional->calcPotentialGGA(df_dt, df_dgamma, grad_t,
                                                     this->derivative, this->max_scale);
            potentialFunction.push_back(pot);
        }
        else {
            FunctionTreeVector<3> df_dgt;
            df_dgt.push_back(xcOutput[2]);
            df_dgt.push_back(xcOutput[3]);
            df_dgt.push_back(xcOutput[4]);
            pot = this->functional->calcPotentialGGA(df_dt, df_dgt,
                                                     this->derivative, this->max_scale);
            potentialFunction.push_back(pot);
        }
    }
    pot = 0;
}

/** @brief clears all data in the XCOperator object
 *
 */
void XCOperator::clear() {
    this->energy = 0.0;
    this->density.clear();
    this->grad_t.clear(true);
    this->grad_a.clear(true);
    this->grad_b.clear(true);
    this->gamma.clear(true);
	this->potentialFunction.clear(true);
    clearApplyPrec();
}

/** @brief given a set of orbitals, it computes the corresponding density function(s)
 *
 * Note: yet another function that does not necessarily belong to this class.
 */
void XCOperator::calcDensity() {
    if (this->orbitals == 0) MSG_ERROR("Orbitals not initialized");
    
    OrbitalVector &phi = *this->orbitals;
    Density &rho = this->density;
    QMPotential &V = *this;
    
    DensityProjector project(this->apply_prec, this->max_scale);
    
    Timer timer1;
    project(rho, phi);
    timer1.stop();
    double t1 = timer1.getWallTime();
    int n1 = rho.getNNodes();
    TelePrompter::printTree(0, "XC density", n1, t1);
    
    if (this->functional->isGGA()) {
        Timer timer2;
        int n2 = calcDensityGradient();
        timer2.stop();
        double t2 = timer2.getWallTime();
        TelePrompter::printTree(0, "XC density gradient", n2, t2);
        printout(1, endl);
    }
    if (this->functional->needsGamma()) calcGamma();
            
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
void XCOperator::calcGamma() {
    FunctionTree<3> * temp;
    if(this->functional->isSpinSeparated()) {
            temp = calcDotProduct(grad_a, grad_a);
            this->gamma.push_back(temp);
            temp = calcDotProduct(grad_a, grad_b);
            this->gamma.push_back(temp);
            temp = calcDotProduct(grad_b, grad_b);
            this->gamma.push_back(temp);
        } else {
            temp = calcDotProduct(grad_t, grad_t);
            this->gamma.push_back(temp);
    }
}

/** @brief computes the gradient of the density
 *
 * For spin-free calculations, the total density is used. For
 * spin-separated calculations both alpha and beta gradients are
 * computed. The results are stored in the correspondig data members
 * of the XCOperator.
 *
 */
int XCOperator::calcDensityGradient() {
    int nNodes = 0;
    if (this->density.isSpinDensity()) {
        grad_a = calcGradient(this->density.alpha());
        grad_b = calcGradient(this->density.beta());
        nNodes += grad_a[0]->getNNodes() + grad_a[1]->getNNodes() + grad_a[2]->getNNodes();
        nNodes += grad_b[0]->getNNodes() + grad_b[1]->getNNodes() + grad_b[2]->getNNodes();
    } else {
        grad_t = calcGradient(this->density.total());
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
FunctionTreeVector<3> XCOperator::calcGradient(FunctionTree<3> &function) {
    if (this->derivative == 0) MSG_ERROR("No derivative operator");
    MWDerivative<3> apply(this->max_scale);
    FunctionTreeVector<3> gradient;
    for (int d = 0; d < 3; d++) {
        FunctionTree<3> *gradient_comp = new FunctionTree<3>(*MRA);
        apply(*gradient_comp, *this->derivative, function, d);
        gradient.push_back(gradient_comp);
    }
    return gradient;
}

/** @brief allocate input arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocared and its pointers assigned to the required input
 * functions.
 *
 */
void XCOperator::setupXCInput() {
   if (this->xcInput != 0) MSG_ERROR("XC input not empty");
    Timer timer;
    println(2, "Preprocessing");

    int nInp = this->functional->getInputLength();
    bool spin = this->functional->isSpinSeparated();
    bool gga = this->functional->isGGA();
    bool gamma = this->functional->needsGamma();

    this->xcInput = allocPtrArray<FunctionTree<3> >(nInp);

    int nUsed = 0;
    nUsed = setupXCInputDensity(nUsed);
    if (gga) {
        nUsed = setupXCInputGradient(nUsed);
    }
    if (nInp != nUsed)  MSG_ERROR("Mismatch between used vs requested");
    for (int i = 0; i < nInp; i++) {
        if (this->xcInput[i] == 0) MSG_ERROR("Invalid XC input");
    }

}

/** @brief sets xcInput pointers for the density
 *
 * Returns the nr. of pointers used for sanity checking
 */
int XCOperator::setupXCInputDensity(int nUsed) {

    bool spinSep = this->functional->isSpinSeparated();
    if(spinSep) {
        this->xcInput[nUsed]     = &this->density.alpha();
        this->xcInput[nUsed + 1] = &this->density.beta();
        nUsed += 2;
    } else {
        this->xcInput[nUsed] = &this->density.total();
        nUsed++;
    }
    return nUsed;
}

/** @brief sets xcInput pointers for the gradient(s)
 *
 * Returns the nr. of pointers used for sanity checking
 */
int XCOperator::setupXCInputGradient(int nUsed) {

    bool spinSep = this->functional->isSpinSeparated();
    bool needsGamma = this->functional->needsGamma();
    
    if(spinSep) {
        if(needsGamma) {
            this->xcInput[nUsed    ] = gamma[0];
            this->xcInput[nUsed + 1] = gamma[1];
            this->xcInput[nUsed + 2] = gamma[2];
            nUsed += 3;
        } else {
            this->xcInput[nUsed    ] = grad_a[0];
            this->xcInput[nUsed + 1] = grad_a[1];
            this->xcInput[nUsed + 2] = grad_a[2];
            this->xcInput[nUsed + 3] = grad_b[0];
            this->xcInput[nUsed + 4] = grad_b[1];
            this->xcInput[nUsed + 5] = grad_b[2];
            nUsed += 6;
        }
    } else {
        if(needsGamma) {
            this->xcInput[nUsed] = gamma[0];
            nUsed += 1;
        } else {
            this->xcInput[nUsed    ] = grad_t[0];
            this->xcInput[nUsed + 1] = grad_t[1];
            this->xcInput[nUsed + 2] = grad_t[2];
            nUsed += 3;
        }
    }
    return nUsed;
}

/** @brief clear the xcInput array
 *
 * the array is just for bookkeeping, therefore it is only necessary
 * to set all pointers to NULL.
 *
 */
void XCOperator::clearXCInput() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");

    int nInp = this->functional->getInputLength();
    bool spin = this->functional->isSpinSeparated();

    for (int i = 0; i < nInp; i++) {
        this->xcInput[i] = 0;
    }
    //LUCA: is this enough? Are we not leaving garbage around?
    this->xcInput = deletePtrArray<FunctionTree<3> >(nInp, &this->xcInput);  
}

/** @brief allocate output arrays for xcfun
 *
 * Based on the xcfun setup, the requested array of FunctionTrees(s)
 * is allocated and the function objects are created, borrowing the
 * grid from the electronic density.
 *
 */
void XCOperator::setupXCOutput() {
    if (this->xcOutput != 0) MSG_ERROR("XC output not empty");
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");
    if (this->xcInput[0] == 0) MSG_ERROR("XC input not initialized");

    GridGenerator<3> grid(this->max_scale);

    // Alloc output trees
    int nOut = this->functional->getOutputLength();
    this->xcOutput = allocPtrArray<FunctionTree<3> >(nOut);

    // Copy grid from input density
    FunctionTree<3> &rho = *this->xcInput[0];
    for (int i = 0; i < nOut; i++) {
        this->xcOutput[i] = new FunctionTree<3>(*MRA);
        grid(*this->xcOutput[i], rho);
    }
}

/** @brief clear the xcOutput array
 *
 * after calling xcfun the array contains intermediate functions which
 * have been employed to obtain the XC potentials. They need to be
 * deleted properly, unless they are used as such (for LDA). In that
 * case the correspondinf pointers are set to NULL when the function
 * "becomes" the potential.
 *
 */
void XCOperator::clearXCOutput() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");

    int nOut = this->functional->getOutputLength();
    this->xcOutput = deletePtrArray<FunctionTree<3> >(nOut, &this->xcOutput);
}

/** @brief evaluation of the functional and its derivatives
 *
 * the data contained in the xcInput is converted in matrix form and fed to the functional.
 * the output matrix is then converted back to function form.
 *
 */
void XCOperator::evaluateXCFunctional() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");
    if (this->xcOutput == 0) MSG_ERROR("XC input not initialized");

    Timer timer;
    println(2, "Evaluating");

    int nInp = this->functional->getInputLength();
    int nOut = this->functional->getOutputLength();

#pragma omp parallel firstprivate(nInp, nOut)
    {
    	int nNodes = this->xcInput[0]->getNEndNodes();
#pragma omp for schedule(guided)
    	for (int n = 0; n < nNodes; n++) {
            MatrixXd inpData, outData;
            compressNodeData(n, nInp, this->xcInput, inpData);
            this->functional->evaluate(this->order, inpData, outData);
            expandNodeData(n, nOut, this->xcOutput, outData);
        }
    }
    for (int i = 0; i < nOut; i++) {
        this->xcOutput[i]->mwTransform(BottomUp);
        this->xcOutput[i]->calcSquareNorm();
    }

    timer.stop();
    double t = timer.getWallTime();
    int n = sumNodes<FunctionTree<3> >(this->xcOutput, nOut);
    TelePrompter::printTree(0, "XC evaluate xcfun", n, t);
    printout(2, endl);
}

/** @brief computes the XC energy as the integral of the functional.
 *
 */
void XCOperator::calcEnergy() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");
    if (this->xcOutput[0] == 0) MSG_ERROR("Invalid XC output");

    Timer timer;
    this->energy = this->xcOutput[0]->integrate();
    timer.stop();
    double t = timer.getWallTime();
    int n = this->xcOutput[0]->getNNodes();
    TelePrompter::printTree(0, "XC energy", n, t);
}

/** @brief scalar product of two FunctionTreeVector(s)
 *
 * param[in] vec_a first vector
 * param[in] vec_b second vector
 *
 * Note: should be a mrcpp functionality.
 *
 */
FunctionTree<3>* XCOperator::calcDotProduct(FunctionTreeVector<3> &vec_a,
                                            FunctionTreeVector<3> &vec_b) {
    if (vec_a.size() != vec_b.size()) MSG_ERROR("Invalid input");

    MWAdder<3> add(-1.0, this->max_scale);
    MWMultiplier<3> mult(-1.0, this->max_scale);
    GridGenerator<3> grid(this->max_scale);

    FunctionTreeVector<3> out_vec;
    for (int d = 0; d < vec_a.size(); d++) {
        FunctionTree<3> &tree_a = vec_a.getFunc(d);
        FunctionTree<3> &tree_b = vec_b.getFunc(d);
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        grid(*out_d, tree_a);
        grid(*out_d, tree_b);
        mult(*out_d, 1.0, tree_a, tree_b, 0);
        out_vec.push_back(out_d);
    }
    FunctionTree<3> *out = new FunctionTree<3>(*MRA);
    grid(*out, out_vec);
    add(*out, out_vec, 0);

    out_vec.clear(true);
    return out;
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
void XCOperator::compressNodeData(int n, int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");
    if (trees[0] == 0) MSG_ERROR("Invalid input");

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
void XCOperator::expandNodeData(int n, int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized output tree " << i);
        VectorXd col_i = data.col(i);
        FunctionNode<3> &node = trees[i]->getEndFuncNode(n);
        node.setValues(col_i);
    }
}

/** @brief fetches the correct index for the potential function to use
 *
 * @param[in] orb the potentialFunction will be applied to this orbital.
 * 
 * Based on the orbital spin, and whether the functional is spin
 * separated, the correct potential index is selected.
 *
 */
int XCOperator::getPotentialFunctionIndex(const Orbital &orb) {
    int orbitalSpin = orb.getSpin();
    bool spinSeparatedFunctional = this->functional->isSpinSeparated();
    int potentialFunctionIndex = -1;
    if (spinSeparatedFunctional and orbitalSpin == Alpha) {
        potentialFunctionIndex = 0;
    }
    else if (spinSeparatedFunctional and orbitalSpin == Beta) {
        potentialFunctionIndex = 1;
    }
    else if (!spinSeparatedFunctional and orbitalSpin == Paired) {
        potentialFunctionIndex = 0;
    }
    else {
        NOT_IMPLEMENTED_ABORT;
    }
    return potentialFunctionIndex;
}
