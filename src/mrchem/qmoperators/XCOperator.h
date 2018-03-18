#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "FunctionTree.h"
#include "FunctionTreeVector.h"
#include "QMPotential.h"
#include "Density.h"

class XCFunctional;
class OrbitalVector;
template<int D> class DerivativeOperator;

/** 
 *  \class XCOperator
 *  \brief Exchange and Correlation operators
 *
 * Ideally, this should handle an arbitrary-order operator (k>1) for
 * response calculations Currently nly the potential is computed
 * (k=1). LDA and GGA functionals are allowed as well as two different ways
 * to compute the XC potentials: either with explicit derivatives or gamma-type derivatives
 * The calss handles the bookkeeping for the input/output of the xcfun library through arrays of
 * FunctionTree(s). 
 *
 *  Testing the output on Sphinx
 *
 *  \author Stig Rune Jensen
 *  \date 2015
 *  
 */
class XCOperator : public QMPotential {
public:
    XCOperator(int k, XCFunctional &F, OrbitalVector &phi, DerivativeOperator<3> *D);
    virtual ~XCOperator();

    double getEnergy() const { return this->energy; }
    void setup(double prec);
    void clear();

    virtual Orbital* operator() (Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;
    
protected:
    const int order;                    ///< Order of kernel derivative
    int nPotentials;                    ///< Number of potential energy functions
    XCFunctional *functional;           ///< External XC functional to be used
    DerivativeOperator<3> *derivative;  ///< External derivative operator
    OrbitalVector *orbitals;            ///< External set of orbitals used to build the density

    double energy;                      ///< XC energy
    Density density;                    ///< Unperturbed density

    FunctionTree<3> **xcInput;          ///< Bookkeeping array to feed XCFun
    FunctionTree<3> **xcOutput;         ///< Bookkeeping array returned by XCFun

    FunctionTreeVector<3> potentialFunction; ///< Storage of the computed potential functions
    FunctionTreeVector<3> grad_a; ///< Gradient of the alpha density        
    FunctionTreeVector<3> grad_b; ///< Gradient of the beta  density        
    FunctionTreeVector<3> grad_t; ///< Gradient of the total density        
    FunctionTreeVector<3> gamma;  ///< Gamma function(s) (three fcns for spin separated calculations)       

    void setupXCInput();
    void setupXCOutput();

    int setupXCInputDensity(int nUsed);
    int setupXCInputGradient(int nUsed);

    void clearXCInput();
    void clearXCOutput();

    void calcDensity();
    int calcDensityGradient();
    void calcGamma();

    void calcPotential();
    bool cropPotential(double prec);

    void calcPotentialLDA();
    
    void calcPotentialGGA();

    void calcEnergy();
    void evaluateXCFunctional();

    void compressTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    void compressNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    FunctionTreeVector<3> calcGradient(FunctionTree<3> &inp);
    FunctionTree<3>* calcDotProduct(FunctionTreeVector<3> &vec_a, FunctionTreeVector<3> &vec_b);

    int getPotentialFunctionIndex(const Orbital & orb);
    
    template<class T>
    int sumNodes(T **trees, int nTrees) {
        int nNodes = 0;
        for (int i = 0; i < nTrees; i++) {
            if (trees[i] != 0) {
                nNodes += trees[i]->getNNodes();
            }
        }
        return nNodes;
    }

    template<class T>
    T** allocPtrArray(int n_funcs) {
        T **ptr = new T*[n_funcs];
        for (int i = 0; i < n_funcs; i++) {
            ptr[i] = 0;
        }
        return ptr;
    }

    template<class T>
    void clearPtrArray(int n_funcs, T **ptr) {
        if (ptr == 0) MSG_FATAL("Clearing NULL pointer");
        for (int i = 0; i < n_funcs; i++) {
            if (ptr[i] != 0) {
                delete ptr[i];
            }
            ptr[i] = 0;
        }
    }

    template<class T>
    T** deletePtrArray(int n_funcs, T ***ptr) {
        if (*ptr != 0) {
            for (int i = 0; i < n_funcs; i++) {
                if ((*ptr)[i] != 0) {
                    delete (*ptr)[i];
                }
                (*ptr)[i] = 0;
            }
            delete[] *ptr;
        }
        return 0;
    }
};


