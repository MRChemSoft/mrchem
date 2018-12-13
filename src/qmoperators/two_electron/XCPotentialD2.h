#pragma once

#include "qmoperators/two_electron/XCPotential.h"

/** @class XCPotential
 *
 * @brief Exchange-Correlation potential defined by a particular (spin) density
 *
 * The XC potential is computed by mapping of the density through a XC functional,
 * provided by the XCFun library. There are two ways of defining the density:
 *
 *  1) Use getDensity() prior to setup() and build the density as you like.
 *  2) Provide a default set of orbitals in the constructor that is used to
 *     compute the density on-the-fly in setup().
 *
 * If a set of orbitals has NOT been given in the constructor, the density
 * MUST be explicitly computed prior to setup(). The density will be computed
 * on-the-fly in setup() ONLY if it is not already available. After setup() the
 * operator will be fixed until clear(), which deletes both the density and the
 * potential.
 *
 * LDA and GGA functionals are supported as well as two different ways to compute
 * the XC potentials: either with explicit derivatives or gamma-type derivatives.
 */

namespace mrchem {

class XCPotentialD2 final : public XCPotential {
public:
    XCPotentialD2(mrdft::XCFunctional *F, OrbitalVector *Phi, OrbitalVector *X, OrbitalVector *Y);
    ~XCPotentialD2();

private:
    OrbitalVector *orbitals_x;               ///< 1st external set of perturbed orbitals used to build the density
    OrbitalVector *orbitals_y;               ///< 2nd external set of perturbed orbitals used to build the density
    mrcpp::FunctionTreeVector<3> potentials; ///< XC Potential functions collected in a vector
    Density *pertDensity_t;                  ///< total first-order perturbed electronic density
    Density *pertDensity_a;                  ///< alpha first-order perturbed electronic density
    Density *pertDensity_b;                  ///< beta  first-order perturbed electronic density

    void setup(double prec) override;
    void clear() override;

    void setupPotential(double prec);
    mrcpp::FunctionTree<3> &getPotential(int orbitalSpin, int densitySpin);
    mrcpp::FunctionTree<3> *buildComponent(int orbital_spin, int density_spin, Density &pert_dens);

    Orbital apply(Orbital phi);

    //LUCA I wanted to include the following declarations in the cpp
    //file but i did not manage to get the syntax (copied from
    //density_utils.copp) right.
    int getPotentialIndex(int orbitalSpin, int densitySpin);
    void setupPerturbedDensity(double prec = -1.0);
};

} //namespace mrchem
