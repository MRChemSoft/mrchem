#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "HelmholtzVector.h"
#include "Orbital.h"
#include "orbital_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param build: precision for construction of Helmholtz operators
 * @param thrs: lambda threshold for re-use of operators
 *
 * The apply_prec is not assigned and no operators are constructed at this point.
 * The setup() function must be called in order to initialize the operators before
 * application.
 */
HelmholtzVector::HelmholtzVector(double build, double thrs)
        : threshold(thrs),
          build_prec(build),
          apply_prec(-1.0) {
}

/** @brief Prepare operators for application
 *
 * @param prec: Apply precision
 * @param Phi: Vector of orbital energies
 *
 * This will construct the necessary HelmholtzOperators from the list of required
 * orbital energies, and assign one operator to each orbital. The orbital energy
 * vector must have the same size as the OrbitalVector on which to apply.
 */
void HelmholtzVector::setup(double prec, const DoubleVector &energies) {
    Printer::printHeader(0, "Setting up Helmholtz operators");
    this->apply_prec = prec;

    Timer timer;
    this->oper_idx.clear();
    this->lambda.clear();
    if (mpi::orb_size > 1) this->clear();
    for (int i = 0; i < energies.size(); i++) {
        double energy = energies(i);
        if (energy > 0.0 ) energy = -0.5;
        int idx = initHelmholtzOperator(energy, i);
        this->oper_idx.push_back(idx);
        if (mpi::orb_size == 1) {
            double mu = (*this)[i].getMu();
            this->lambda.push_back(-0.5*mu*mu);
        } else {
            this->lambda.push_back(energy);
        }
    }
    if (mpi::orb_size == 1) clearUnused();

    timer.stop();
    Printer::printFooter(0, timer, 2);
}

/** @brief Initialize a Helmholtz operator with given energy
 *
 * @param energy: Orbital energy
 * @param i: Orbital index (position in OrbitalVector)
 *
 * This will construct a singel HelmholtzOperator with the given energy, and
 * assign the operator to a specific orbital. If re-use is activated we will
 * first look through the existing operators and assign to one of them if
 * the energy matches.
 */
int HelmholtzVector::initHelmholtzOperator(double energy, int i) {
    if (energy > 0.0) MSG_ERROR("Complex Helmholtz not available");
    double mu = std::sqrt(-2.0*energy);
    if (mpi::orb_size == 1) {
        for (int j = 0; j < this->operators.size(); j++) {
            double mu_j = this->operators[j]->getMu();
            double muDiff = mu - mu_j;
            double relDiff = muDiff/mu;
            if (std::abs(relDiff) < this->threshold and false) {
                double l = -0.5*mu_j*mu_j;
                Printer::printDouble(0, "Re-using operator with lambda", l, 5);
                return j;
            }
        }
    }
    mrcpp::HelmholtzOperator *oper = 0;
    if( i%mpi::orb_size == mpi::orb_rank) {
        Printer::printDouble(0, "Creating operator with lambda", energy, 5);
        oper = new mrcpp::HelmholtzOperator(*MRA, mu, this->build_prec);
    }
    this->operators.push_back(oper);

    return this->operators.size() - 1;
}

/** @brief Clear operators after application
 *
 * This will clear all existing operators, even if re-use is activated.
 */
void HelmholtzVector::clear() {
    int nOper = this->operators.size();
    for (int j = 0; j < nOper; j++) {
        if (this->operators[j] != 0) {
            delete this->operators[j];
            this->operators[j] = 0;
        }
    }
    this->operators.clear();
    this->oper_idx.clear();
    this->lambda.clear();
}

/** @brief Clear unused operators
 *
 * This will clear all operators that currently doesn't point to any orbital.
 */
void HelmholtzVector::clearUnused() {
    std::vector<mrcpp::HelmholtzOperator *> tmp;
    int nIdx = this->oper_idx.size();
    int nOper = this->operators.size();
    IntVector nUsed = IntVector::Zero(nOper);
    for (int i = 0; i < nIdx; i++) {
        int n = this->oper_idx[i];
        nUsed(n) = nUsed(n) + 1;
    }
    for (int i = 0; i < nOper; i++) {
        if (nUsed(i) == 0) {
            delete this->operators[i];
        } else {
            tmp.push_back(this->operators[i]);
            int newIdx = tmp.size() - 1;
            for (int j = 0; j < nIdx; j++) {
                if (this->oper_idx[j] == i) this->oper_idx[j] = newIdx;
            }
        }
        this->operators[i] = 0;
    }
    this->operators.clear();
    for (int i = 0; i < tmp.size(); i++) {
        this->operators.push_back(tmp[i]);
        tmp[i] = 0;
    }
}

/** @brief Return the operator corresponding to the i-th orbital */
mrcpp::HelmholtzOperator& HelmholtzVector::operator[](int i) {
    int idx = this->oper_idx[i];
    if (idx < 0 or idx >= operators.size()) MSG_ERROR("Invalid operator index");
    mrcpp::HelmholtzOperator *oper = this->operators[idx];
    if (oper == 0) MSG_ERROR("Operator null pointer");
    return *oper;
}

/** @brief Return the operator corresponding to the i-th orbital */
const mrcpp::HelmholtzOperator& HelmholtzVector::operator[](int i) const {
    int idx = this->oper_idx[i];
    if (idx < 0 or idx >= operators.size()) MSG_ERROR("Invalid operator index");
    const mrcpp::HelmholtzOperator *oper = this->operators[idx];
    if (oper == 0) MSG_ERROR("Operator null pointer");
    return *oper;
}

/** @brief Return a vector of the lambda parameters of the currently operators*/
DoubleVector HelmholtzVector::getLambdaVector() const {
    int nLambda = this->lambda.size();
    DoubleVector L = DoubleVector::Zero(nLambda);
    for (int i = 0; i < nLambda; i++) {
        L(i) = this->lambda[i];
    }
    return L;
}

/** @brief Return a diagonal matrix of the lambda parameters of the currently operators*/
ComplexMatrix HelmholtzVector::getLambdaMatrix() const {
    ComplexVector lambda = getLambdaVector().cast<ComplexDouble>();
    return lambda.asDiagonal();
}

/** @brief Prints the total number of trees and nodes kept in the operator set */
int HelmholtzVector::printTreeSizes() const {
    int totNodes = 0;
    int totTrees = 0;
    int nOperators = this->operators.size();
    for (int i = 0; i < nOperators; i++) {
        if(this->operators[i] != 0){
            int nTrees = this->operators[i]->size();
            for (int j = 0; j < nTrees; j++) {
                totNodes += this->operators[i]->getComponent(j).getNNodes();
            }
            totTrees += nTrees;
        }
    }
    return totNodes;
}

/** @brief Apply the i-th operator to an orbital
 *
 * @param i: Orbital index (position in OrbtialVector)
 * @param inp: Orbital on which to apply
 *
 * This will apply the HelmholtzOperator onto both the real and imaginary
 * parts of the input orbital.
 */
Orbital HelmholtzVector::operator()(int i, Orbital inp) {
    mrcpp::HelmholtzOperator &H_i = (*this)[i];

    Orbital out = inp.paramCopy();
    if (inp.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(this->apply_prec, out.real(), H_i, inp.real());
    }
    if (inp.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(this->apply_prec, out.imag(), H_i, inp.imag());
        if (inp.conjugate()) out.imag().rescale(-1.0);
    }
    return out;
}

/** @brief Apply all operators to the corresponding orbital
 *
 * @param inp: Orbitals on which to apply
 *
 * This will produce an output OrbitalVector where each of the HelmholtzOperators
 * have been applied to the corresponding orbital in the input vector.
 */
OrbitalVector HelmholtzVector::operator()(OrbitalVector &inp) {
    Printer::printHeader(0, "Applying Helmholtz operators");
    println(0, " Orb    RealNorm   Nodes     ImagNorm   Nodes     Timing");
    Printer::printSeparator(0, '-');
    int oldprec = Printer::setPrecision(5);

    OrbitalVector out = orbital::param_copy(inp);
    HelmholtzVector &H = (*this);

    Timer tottimer;
    for (int i = 0; i < inp.size(); i++) {
        if (not mpi::my_orb(inp[i])) continue;

        Timer timer;
        out[i] = H(i, inp[i]);

        int rNodes = out[i].getNNodes(NUMBER::Real);
        int iNodes = out[i].getNNodes(NUMBER::Imag);
        double rNorm = 0.0;
        double iNorm = 0.0;
        if (out[i].hasReal()) rNorm = std::sqrt(out[i].real().getSquareNorm());
        if (out[i].hasImag()) iNorm = std::sqrt(out[i].imag().getSquareNorm());

        timer.stop();
        Printer::setPrecision(5);
        printout(0, std::setw(3) << i);
        printout(0, " " << std::setw(14) << rNorm);
        printout(0, " " << std::setw(5) << rNodes);
        printout(0, " " << std::setw(14) << iNorm);
        printout(0, " " << std::setw(5) << iNodes);
        printout(0, std::setw(14) << timer.getWallTime() << std::endl);
    }

    mpi::barrier(mpi::comm_orb); // barrier to align printing

    tottimer.stop();
    Printer::printFooter(0, tottimer, 2);
    Printer::setPrecision(oldprec);
    return out;
}

} //namespace mrchem
