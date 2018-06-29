#pragma once

#include "qmoperators.h"
#include "qmfunction_utils.h"

/** @class HelmholtzVector
 *
 * @brief Container of HelmholtzOperators for a corresponding OrbtialVector
 *
 * This class assigns one HelmholtzOperator to each orbital in an OrbitalVector.
 * The operators can be re-used if several orbitals share the same (or similar)
 * energy, or if the change in energy is small relative to the previous iteration.
 */

namespace mrchem {

class HelmholtzVector final {
public:
    HelmholtzVector(double build, double thrs = -1.0);

    void setup(double prec, const DoubleVector &energies);
    void clear();

    void setThreshold(double thrs) { this->threshold = thrs; }
    double getThreshold() const { return this->threshold; }

    double getLambda(int i) const { return this->lambda[i]; }
    DoubleVector getLambdaVector() const;
    ComplexMatrix getLambdaMatrix() const;

    mrcpp::HelmholtzOperator& operator[](int i);
    const mrcpp::HelmholtzOperator& operator[](int i) const;

    int printTreeSizes() const;

    Orbital operator()(int i, Orbital inp);
    OrbitalVector operator()(OrbitalVector &inp);

private:
    double threshold;   ///< For re-using operators. Negative means always recreate
    double build_prec;  ///< Precision for construction of Helmholtz operators
    double apply_prec;  ///< Precision for application of Helmholtz operators

    std::vector<int> oper_idx;  ///< Points to a HelmholtzOperator in the operators vector
    std::vector<double> lambda; ///< The lambda value used for the corresponding HelmholtzOperator
    std::vector<mrcpp::HelmholtzOperator *> operators; ///< Vector of Helmholtz operators

    int initHelmholtzOperator(double energy, int i);
    void clearUnused();
};

} //namespace mrchem
