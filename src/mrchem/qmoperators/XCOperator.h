#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "QMPotential.h"
#include "Density.h"

class XCFunctional;
class OrbitalVector;
template<int D> class FunctionTree;
template<int D> class FunctionTreeVector;
template<int D> class DerivativeOperator;

/** 
 *  \class XCOperator
 *  \brief Interface class to compute DFT functionaals and derivatives
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

protected:
    const int order;                    ///< Order of kernel derivative
    XCFunctional *functional;           ///< Pointer to external object
    DerivativeOperator<3> *derivative;  ///< Pointer to external object
    OrbitalVector *orbitals;            ///< Pointer to external object

    double energy;                      ///< XC energy
    Density density;                    ///< Unperturbed density
    Density gradient[3];                ///< Unperturbed density gradient

    FunctionTree<3> **xcInput;          ///< XCFun input
    FunctionTree<3> **xcOutput;         ///< XCFun output

    void setupXCInput();
    void setupXCOutput();

    void clearXCInput();
    void clearXCOutput();

    void calcDensity();
    void calcDensityGradient(Density *dRho, Density &rho);

    virtual void calcPotential() = 0;
    bool cropPotential(double prec);

    void calcEnergy();
    void evaluateXCFunctional();

    void compressTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandTreeData(int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    void compressNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, FunctionTree<3> **trees, Eigen::MatrixXd &data);

    FunctionTreeVector<3> calcGradient(FunctionTree<3> &inp);
    FunctionTree<3>* calcDivergence(FunctionTreeVector<3> &inp);
    FunctionTree<3>* calcDotProduct(FunctionTreeVector<3> &vec_a, FunctionTreeVector<3> &vec_b);
    FunctionTree<3>* calcGradDotPotDensVec(FunctionTree<3> &V, FunctionTreeVector<3> &rho);

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


