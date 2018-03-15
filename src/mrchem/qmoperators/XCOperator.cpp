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
    std::cout << "Is spin sep 2" << spin << std::endl;
    nPotentials = spin ? k + 1 : 1; /// k+1 potentials if spin separated, otherwise just one.
    density.setIsSpinDensity(spin);
    gradient[0].setIsSpinDensity(spin);
    gradient[1].setIsSpinDensity(spin);
    gradient[2].setIsSpinDensity(spin);
    F.evalSetup(k);
}

XCOperator::~XCOperator() {
    if (this->xcInput != 0) MSG_ERROR("XC input not deallocated");
    if (this->xcOutput != 0) MSG_ERROR("XC output not deallocated");
    this->functional = 0;
    this->derivative = 0;
}

/** \brief Provides the right sequence of internal operations to compute the XC potential 
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

Orbital* XCOperator::operator() (Orbital &phi) {
    FunctionTree<3> * potential = this->potentialFunction[this->getPotentialFunctionIndex(phi)];
    this->setReal(potential);
    this->setImag(0);
    Orbital * Vphi = QMPotential::operator()(phi); 
    this->clearReal(false);
    this->clearImag(false);
    return Vphi;
}

Orbital* XCOperator::adjoint(Orbital &phi) {
    NOT_IMPLEMENTED_ABORT;
}

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

void XCOperator::calcPotentialLDA() {

    if (this->order != 1) {
        NOT_IMPLEMENTED_ABORT;
    }
    for (int i = 0; in < this->nPotentials) {
        int outputIndex = i + 1;
        if (xcOutput[outputIndex] == 0) MSG_ERROR("Invalid XC output");
        potentialFunction.push_back(xcOutput[outputIndex]);
        xcOutput[outputIndex] = 0;
    }
}


void XCOperator::calcPotentialGGA() {
    
    bool spin = this->functional->isSpinSeparated();
    bool gamma = this->functional->needsGamma();
    if(spin) {
        FunctionTree & df_da   = xcOutput[1];
        FunctionTree & df_db   = xcOutput[2];
        if(gamma) {
            FunctionTree & df_dgaa = xcOutput[3];
            FunctionTree & df_dgab = xcOutput[4];
            FunctionTree & df_dgbb = xcOutput[5];
            pot = calcPotentialGGA(df_da, df_db, df_dgaa, df_dgab, grad_a, grad_b);
            potentialFunction.push_back(pot);
            pot = calcPotentialGGA(df_db, df_da, df_dgbb, df_dgab, grad_b, grad_a);
            potentialFunction.push_back(pot);
        }
        else {
            FunctionTreeVector<3> df_dga;
            FunctionTreeVector<3> df_dgb;
            df_dga.push_back(xcOutput[3]);
            df_dga.push_back(xcOutput[4]);
            df_dga.push_back(xcOutput[5]);
            df_dga.push_back(xcOutput[6]);
            df_dga.push_back(xcOutput[7]);
            df_dga.push_back(xcOutput[8]);
            pot = calcPotentialGGA(df_da, df_dga);
            potentialFunction.push_back(pot);
            pot = calcPotentialGGA(df_db, df_dgb);
            potentialFunction.push_back(pot);
        }
            
    }
    else {
        FunctionTree & df_dt   = xcOutput[1];
        if(gamma) {
            FunctionTree & df_dgamma   = xcOutput[2];
            pot = calcPotentialGGA(df_dt, df_dgamma, grad_t);
            potentialFunction.push_back(pot);
        }
        else {
            FunctionTreeVector<3> df_dgt;
            df_dga.push_back(xcOutput[3]);
            df_dga.push_back(xcOutput[4]);
            df_dga.push_back(xcOutput[5]);
            pot = calcPotentialGGA(df_dt, df_dgt);
            potentialFunction.push_back(pot);
        }
    }
}

/** \brief Cleanup function to call after the the XC potential has been calculated
 */
void XCOperator::clear() {
    this->energy = 0.0;
    this->density.clear();
    this->gradient[0].clear();
    this->gradient[1].clear();
    this->gradient[2].clear();
	for (int i = 0; i < nPotentials; i++) {
		this->potentialFunction[i]->clear();
		this->potentialFunction.pop_back();
	}
    clearApplyPrec();
}

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
        n2 = calcDensityGradient(rho);
        timer2.stop();
        double t2 = timer2.getWallTime();
        TelePrompter::printTree(0, "XC density gradient", n2, t2);
        printout(1, endl);
    }
}

int XCOperator::calcDensityGradient(Density &rho) {
    if (rho.isSpinDensity()) {
        grad_a = calcGradient(rho.alpha());
        grad_b = calcGradient(rho.beta());
        nNodes  = grad_a[0]->getNNodes() + grad_a[1]->getNNodes() + grad_a[2]->getNNodes()
        nNodes += grad_b[0]->getNNodes() + grad_b[1]->getNNodes() + grad_b[2]->getNNodes()
    } else {
        grad_t = calcGradient(rho.total());
        nNodes = grad_t[0]->getNNodes() + grad_t[1]->getNNodes() + grad_t[2]->getNNodes()
    }
    return nNodes;
}

FunctionTreeVector<3> XCOperator::calcGradient(FunctionTree<3> &inp) {
    if (this->derivative == 0) MSG_ERROR("No derivative operator");
    MWDerivative<3> apply(this->max_scale);
    
    FunctionTreeVector<3> out;
    for (int d = 0; d < 3; d++) {
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        apply(*out_d, *this->derivative, inp, d);
        out.push_back(out_d);
    }
    return out;
}

void XCOperator::setupXCInput() {
   if (this->xcInput != 0) MSG_ERROR("XC input not empty");
    Timer timer;
    println(2, "Preprocessing");

    int nInp = this->functional->getInputLength();
    bool spin = this->functional->isSpinSeparated();
    bool gga = this->functional->isGGA();
    bool gamma = this->functional->needsGamma();

    std::cout << "Input length " << nInp << std::endl;
    
    this->xcInput = allocPtrArray<FunctionTree<3> >(nInp);

    int nUsed = 0;
    nUsed = setupXCInputDensity(nUsed, spin);
    if (gga) {
        nUsed = setupXCInputGradient(nUsed, spin, gamma);
    }
    
    // sanity check
    if (nInp != nUsed)  MSG_ERROR("Mismatch between used vs requested");
    for (int i = 0; i < nInp; i++) {
        if (this->xcInput[i] == 0) MSG_ERROR("Invalid XC input");
    }

}

int XCOperator::setupXCInputDensity(int nUsed, bool spin) {
    if(spin) {
        this->xcInput[nUsed]     = &this->density.alpha();
        this->xcInput[nUsed + 1] = &this->density.beta();
        nUsed += 2;
    } else {
        this->xcInput[nUsed] = &this->density.total();
        nUsed++;
    }
    return nUsed;
}

int XCOperator::setupXCInputGradient(int nUsed, bool spin, bool gamma) {

    FunctionTreeVector<3> vec_a, vec_b, vec_t;
    
    vec_a.push_back(&gradient[0].alpha());
    vec_a.push_back(&gradient[1].alpha());
    vec_a.push_back(&gradient[2].alpha());

    vec_b.push_back(&gradient[0].beta());
    vec_b.push_back(&gradient[1].beta());
    vec_b.push_back(&gradient[2].beta());
                    
    vec_t.push_back(&gradient[0].total());
    vec_t.push_back(&gradient[1].total());
    vec_t.push_back(&gradient[2].total());
    
    if(spin) {
        if(gamma) {
            this->xcInput[nUsed    ] = calcDotProduct(vec_a, vec_a);
            this->xcInput[nUsed + 1] = calcDotProduct(vec_a, vec_b);
            this->xcInput[nUsed + 2] = calcDotProduct(vec_b, vec_b);
            nUsed += 3;
        } else {
            this->xcInput[nUsed    ] = vec_a[0];
            this->xcInput[nUsed + 1] = vec_a[1];
            this->xcInput[nUsed + 2] = vec_a[2];
            this->xcInput[nUsed + 3] = vec_b[0];
            this->xcInput[nUsed + 4] = vec_b[1];
            this->xcInput[nUsed + 5] = vec_b[2];
            nUsed += 6;
        }
    } else {
        if(gamma) {
            this->xcInput[nUsed] = calcDotProduct(vec_t, vec_t);
            nUsed += 1;
        } else {
            this->xcInput[nUsed    ] = vec_t[0];
            this->xcInput[nUsed + 1] = vec_t[1];
            this->xcInput[nUsed + 2] = vec_t[2];
            nUsed += 3;
        }
    }
    vec_a.clear();
    vec_b.clear();
    vec_t.clear();
    
    return nUsed;
}

void XCOperator::clearXCInput() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");

    int nInp = this->functional->getInputLength();
    bool spin = this->functional->isSpinSeparated();

    // these belong to density
    this->xcInput[0] = 0;
    if (spin) this->xcInput[1] = 0;

    // the rest should be deleted
    this->xcInput = deletePtrArray<FunctionTree<3> >(nInp, &this->xcInput);  //LUCA: is this enough? Are we not leaving garbage around?
}

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

void XCOperator::clearXCOutput() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");

    int nOut = this->functional->getOutputLength();
    this->xcOutput = deletePtrArray<FunctionTree<3> >(nOut, &this->xcOutput);
}

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

void XCOperator::calcEnergy() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");
    if (this->xcOutput[0] == 0) MSG_ERROR("Invalid XC output");

    Timer timer;
    this->energy = this->xcOutput[0]->integrate();
    timer.stop();
    double t = timer.getWallTime();
    int n = this->xcOutput[0]->getNNodes();
    std::cout << "XC energy " << this->energy << std::endl;
    TelePrompter::printTree(0, "XC energy", n, t);
}

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

void XCOperator::compressTreeData(int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");
    if (trees[0] == 0) MSG_ERROR("Invalid input");

    FunctionTree<3> &tree = *trees[0];
    int nCoefs = tree.getTDim()*tree.getKp1_d()*tree.getNEndNodes();
    data = MatrixXd::Zero(nCoefs, nFuncs);

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized input tree");
        VectorXd col_i;
        trees[i]->getEndValues(col_i);
        data.col(i) = col_i;
    }
}

void XCOperator::expandTreeData(int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized output tree " << i);
        VectorXd col_i = data.col(i);
        trees[i]->setEndValues(col_i);
    }
}

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

void XCOperator::expandNodeData(int n, int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized output tree " << i);
        VectorXd col_i = data.col(i);
        FunctionNode<3> &node = trees[i]->getEndFuncNode(n);
        node.setValues(col_i);
    }
}

int XCOperator::getPotentialFunctionIndex(const Orbital &orb) {
    int orbitalSpin = orb.getSpin();
    bool spinSeparatedFunctional = this->functional->isSpinSeparated();
    int orbitalOccupancy = orb.getOccupancy();
    int potentialFunctionIndex = -1;
    std::cout << "orbSpin " << orbitalSpin
              << "ssf "     << spinSeparatedFunctional
              << "oo "      << orbitalOccupancy << std::endl; 
    if (spinSeparatedFunctional && orbitalSpin == Alpha) {
        potentialFunctionIndex = 0;
    }
    else if (spinSeparatedFunctional && orbitalSpin == Beta) {
        potentialFunctionIndex = 1;
    }
    else if (!spinSeparatedFunctional && orbitalSpin == Paired) {
        potentialFunctionIndex = 0;
    }
    else {
        NOT_IMPLEMENTED_ABORT;
    }
    return potentialFunctionIndex;
}
