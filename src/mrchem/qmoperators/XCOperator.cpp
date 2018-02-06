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

void XCOperator::calcDensity() {
    if (this->orbitals == 0) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &phi = *this->orbitals;
    Density &rho = this->density;
    Density *dRho = &this->gradient[0];
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
        calcDensityGradient(dRho, rho);
        timer2.stop();
        double t2 = timer2.getWallTime();
        int n2 = 0;
        n2 += dRho[0].getNNodes();
        n2 += dRho[1].getNNodes();
        n2 += dRho[2].getNNodes();
        TelePrompter::printTree(0, "XC density gradient", n2, t2);
        printout(1, endl);
    }
}

void XCOperator::calcDensityGradient(Density *dRho, Density &rho) {
    if (rho.isSpinDensity()) {
        FunctionTreeVector<3> grad_a = calcGradient(rho.alpha());
        dRho[0].setDensity(Density::Alpha, grad_a[0]);
        dRho[1].setDensity(Density::Alpha, grad_a[1]);
        dRho[2].setDensity(Density::Alpha, grad_a[2]);

        FunctionTreeVector<3> grad_b = calcGradient(rho.beta());
        dRho[0].setDensity(Density::Beta, grad_b[0]);
        dRho[1].setDensity(Density::Beta, grad_b[1]);
        dRho[2].setDensity(Density::Beta, grad_b[2]);
    } else {
        FunctionTreeVector<3> grad_t = calcGradient(rho.total());
        dRho[0].setDensity(Density::Total, grad_t[0]);
        dRho[1].setDensity(Density::Total, grad_t[1]);
        dRho[2].setDensity(Density::Total, grad_t[2]);
    }
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

FunctionTree<3>* XCOperator::calcDivergence(FunctionTreeVector<3> &inp) {
    if (this->derivative == 0) MSG_ERROR("No derivative operator");
    MWAdder<3> add(-1.0, this->max_scale);
    MWDerivative<3> apply(this->max_scale);
    GridGenerator<3> grid(this->max_scale);

    FunctionTreeVector<3> tmp_vec;
    for (int d = 0; d < 3; d++) {
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        apply(*out_d, *this->derivative, *inp[d], d);
        tmp_vec.push_back(out_d);
    }
    FunctionTree<3> *out = new FunctionTree<3>(*MRA);
    grid(*out, tmp_vec);
    add(*out, tmp_vec, 0); // Addition on union grid
    tmp_vec.clear(true);
    return out;
}

/** Compute the required input functions for XCFun. In the case of GGA
 * the spin density gradients are computed. Then densities and gradients
 * are sorted in the input function array in the correct order for XCFun.
 * We define \f$ \gamma_{\alpha\beta} = \nabla\rho_\alpha\cdot\nabla\rho_\beta\f$
 * LDA: \f$ (\rho) \f$
 * GGA: \f$ (\rho, \gamma) \f$
 * Spin LDA: \f$ (\rho_\alpha,\rho_\beta) \f$
 * Spin GGA: \f$ (\rho_\alpha,\rho_\beta, \gamma_{\alpha\alpha},
 * \gamma_{\alpha\beta}, \gamma_{\beta\beta}) \f$
 */
void XCOperator::setupXCInput() {
    if (this->xcInput != 0) MSG_ERROR("XC input not empty");
    Timer timer;
    println(2, "Preprocessing");

    int nInp = this->functional->getInputLength();
    bool spin = this->functional->isSpinSeparated();
    bool gga = this->functional->isGGA();

    Density &rho = this->density;
    Density &rho_x = this->gradient[0];
    Density &rho_y = this->gradient[1];
    Density &rho_z = this->gradient[2];

    this->xcInput = allocPtrArray<FunctionTree<3> >(nInp);

    if (not spin) {
        this->xcInput[0] = &rho.total();
        if (gga) {
            FunctionTreeVector<3> vec;
            vec.push_back(&rho_x.total());
            vec.push_back(&rho_y.total());
            vec.push_back(&rho_z.total());
            this->xcInput[1] = calcDotProduct(vec, vec);
            vec.clear();
        }
    } else {
        this->xcInput[0] = &rho.alpha();
        this->xcInput[1] = &rho.beta();
        if (gga) {
            FunctionTreeVector<3> vec;
            vec.push_back(&rho_x.alpha());
            vec.push_back(&rho_y.alpha());
            vec.push_back(&rho_z.alpha());
            this->xcInput[1] = calcDotProduct(vec, vec);
            vec.clear();

            vec.push_back(&rho_x.beta());
            vec.push_back(&rho_y.beta());
            vec.push_back(&rho_z.beta());
            this->xcInput[2] = calcDotProduct(vec, vec);
            vec.clear();
        }
    }

    // sanity check
    for (int i = 0; i < nInp; i++) {
        if (this->xcInput[i] == 0) MSG_ERROR("Invalid XC input");
    }

    timer.stop();
    double t = timer.getWallTime();
    int n = sumNodes<FunctionTree<3> >(this->xcInput, nInp);
    TelePrompter::printTree(0, "XC preprocess xcfun", n, t);
    printout(2, endl);
}

void XCOperator::clearXCInput() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");

    int nInp = this->functional->getInputLength();
    bool spin = this->functional->isSpinSeparated();

    // these belong to density
    this->xcInput[0] = 0;
    if (spin) this->xcInput[1] = 0;

    // the rest should be deleted
    this->xcInput = deletePtrArray<FunctionTree<3> >(nInp, &this->xcInput);
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
    TelePrompter::printTree(0, "XC energy", n, t);
}

FunctionTree<3>* XCOperator::calcGradDotPotDensVec(FunctionTree<3> &V, FunctionTreeVector<3> &rho) {
    MWMultiplier<3> mult(-1.0, this->max_scale);
    GridGenerator<3> grid(this->max_scale);

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
    FunctionTree<3> *result = calcDivergence(vec);
    vec.clear(true);

    timer.stop();
    double t = timer.getWallTime();
    int n = result->getNNodes();
    TelePrompter::printTree(2, "Gradient", n, t);
    return result;
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

/*
Orbital* XCOperator::operator() (Orbital &orb_p) {
    QMPotential &V_p = this->potential[0];
    QMPotential &V_a = this->potential[1];
    QMPotential &V_b = this->potential[2];

    if (orb_p.getSpin() == Density::Total) return V_p(orb_p);
    if (orb_p.getSpin() == Density::Alpha) return V_a(orb_p);
    if (orb_p.getSpin() == Density::Beta) return V_b(orb_p);

    MSG_ERROR("Invalid spin");
    return 0;
}
*/

/*
Orbital* XCOperator::adjoint(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
}
*/

/*
int XCOperator::printTreeSizes() const {
    NOT_IMPLEMENTED_ABORT;
    int nNodes = 0;
    nNodes += this->density_0.printTreeSizes();
    this->gradient_0[0]->printTreeSizes();
    this->gradient_0[1]->printTreeSizes();
    this->gradient_0[2]->printTreeSizes();

    if (this->potential != 0) {
        nNodes += this->potential[0]->printTreeSizes();
        nNodes += this->potential[1]->printTreeSizes();
        nNodes += this->potential[2]->printTreeSizes();
    }
    int nInput = 0;
    int inTrees = 0;
    if (this->xcInput != 0) {
        int nFuncs = this->functional->getInputLength();
        for (int i = 0; i < nFuncs; i++) {
            if (this->xcInput[i] != 0) {
                nInput += this->xcInput[i]->getNNodes();
                inTrees++;
            }
        }
    }
    int nOutput = 0;
    int outTrees = 0;
    if (this->xcOutput != 0) {
        int nFuncs = this->functional->getOutputLength();
        for (int i = 0; i < nFuncs; i++) {
            if (this->xcOutput[i] != 0) {
                nOutput += this->xcOutput[i]->getNNodes();
                outTrees++;
            }
        }
    }
    println(0, " XC input          " << setw(15) << inTrees << setw(25) << nInput);
    println(0, " XC output         " << setw(15) << outTrees << setw(25) << nOutput);

    nNodes += (nInput + nOutput);
    return nNodes;
}
*/

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

