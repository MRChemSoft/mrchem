#pragma once

#include <Eigen/Core>

#include "MRCPP/MWFunctions"
#include "XCFun/xcfun.h"

/** 
 *  @class XCFunctional
 *  @brief Compute XC functional with XCFun
 *
 *  Interface class for the XCFun library
 *
 * Depending on the mode chosen, xcfun needs either the gamma
 * functions or the explicit gradients. The first mode is possibly
 * more efficient (fewer functions to compute/handle), whereas the
 * other is simpler to implement. We keep both options open and
 * compute the gradient invariants if and when necessary.
 *
 *
 */

namespace mrdft {

enum class DensityType { Total, Alpha, Beta };

class XCFunctional final {
public:
    XCFunctional(mrcpp::MultiResolutionAnalysis<3> &mra, bool spin);
    ~XCFunctional();

    bool hasDensity() const;
    mrcpp::FunctionTree<3> & getDensity(DensityType type);

    void setDensityCutoff(double thrs) { this->cutoff = thrs; }
    void setFunctional(const std::string &name, double coef = 1.0);

    void setUseGamma(bool use) { this->use_gamma = use; }
    bool useGamma() const { return this->use_gamma; }

    bool isLDA() const { return (not (isGGA() || isMetaGGA())); }
    bool isGGA() const { return (xc_is_gga(this->functional)); }
    bool isMetaGGA() const { return (xc_is_metagga(this->functional)); }
    bool isSpinSeparated() const { return this->spin_separated; }

    void evalSetup(int order);
    void setup();
    void clear();

    void evaluate();
    double calcEnergy();
    mrcpp::FunctionTreeVector<3> calcPotential();

 protected:
    const bool spin_separated;                  ///< Spin polarization
    const mrcpp::MultiResolutionAnalysis<3> MRA;///< Computational domain

    bool use_gamma;                             ///< Whether gamma-type or explicit derivatives are used
    double cutoff;                              ///< Below the cutoff value, the density will be considered zero

    xc_functional functional;                   ///< The functional in the XCFun library (struct from xcfun library)
    mrcpp::DerivativeOperator<3> *derivative;   ///< Derivative operator

    mrcpp::FunctionTree<3> *rho_a;              ///< Alpha density
    mrcpp::FunctionTree<3> *rho_b;              ///< Beta density
    mrcpp::FunctionTree<3> *rho_t;              ///< Total density
    mrcpp::FunctionTreeVector<3> grad_a;        ///< Gradient of the alpha density
    mrcpp::FunctionTreeVector<3> grad_b;        ///< Gradient of the beta  density
    mrcpp::FunctionTreeVector<3> grad_t;        ///< Gradient of the total density
    mrcpp::FunctionTreeVector<3> gamma;         ///< Gamma function(s) (three fcns for spin separated calculations)

    mrcpp::FunctionTreeVector<3> xcInput;       ///< Bookkeeping array to feed XCFun
    mrcpp::FunctionTreeVector<3> xcOutput;      ///< Bookkeeping array returned by XCFun

    int getInputLength() const { return xc_input_length(this->functional); }
    int getOutputLength() const { return xc_output_length(this->functional); }

    void setupXCInput();
    void setupXCOutput();
    int setupXCInputDensity();
    int setupXCInputGradient();

    void evaluateBlock(Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const;
    void compressNodeData(int n, int nFuncs, mrcpp::FunctionTreeVector<3> trees, Eigen::MatrixXd &data);
    void expandNodeData(int n, int nFuncs, mrcpp::FunctionTreeVector<3> trees, Eigen::MatrixXd &data);

    void calcPotentialLDA(mrcpp::FunctionTreeVector<3> &potentials);
    void calcPotentialGGA(mrcpp::FunctionTreeVector<3> &potentials);

    mrcpp::FunctionTree<3> * calcPotentialGGA(mrcpp::FunctionTree<3> & df_drho,
                                              mrcpp::FunctionTreeVector<3> & df_dgr);
    mrcpp::FunctionTree<3> * calcPotentialGGA(mrcpp::FunctionTree<3> & df_drho,
                                              mrcpp::FunctionTree<3> & df_dgamma,
                                              mrcpp::FunctionTreeVector<3> grad_rho);
    mrcpp::FunctionTree<3> * calcPotentialGGA(mrcpp::FunctionTree<3> & df_drhoa,
                                              mrcpp::FunctionTree<3> & df_dgaa,
                                              mrcpp::FunctionTree<3> & df_dgab,
                                              mrcpp::FunctionTreeVector<3> grad_rhoa,
                                              mrcpp::FunctionTreeVector<3> grad_rhob);

    mrcpp::FunctionTree<3> * addPotentialContributions(mrcpp::FunctionTreeVector<3> & contributions);
    mrcpp::FunctionTree<3> * calcGradDotPotDensVec(mrcpp::FunctionTree<3> &V, mrcpp::FunctionTreeVector<3> &rho);
    mrcpp::FunctionTree<3> * calcDotProduct(mrcpp::FunctionTreeVector<3> &vec_a, mrcpp::FunctionTreeVector<3> &vec_b);
    mrcpp::FunctionTree<3> * calcDivergence(mrcpp::FunctionTreeVector<3> &inp);
    mrcpp::FunctionTreeVector<3> calcGradient(mrcpp::FunctionTree<3> &inp);
};
 
} //namespace mrdft
