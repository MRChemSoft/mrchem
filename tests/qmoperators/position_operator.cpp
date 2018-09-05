#include "catch.hpp"

#include "mrchem.h"

#include "HarmonicOscillatorFunction.h"
#include "PositionOperator.h"
#include "Orbital.h"
#include "orbital_utils.h"

using namespace mrchem;
using namespace orbital;

namespace position_operator {

TEST_CASE("PositionOperator", "[position_operator]") {
    const double prec = 1.0e-3;
    const double thrs = prec*prec;

    int nFuncs = 3;
    OrbitalVector Phi;
    for (int n = 0; n < nFuncs; n++) {
        int nu[3] = {n,0,0};
        HarmonicOscillatorFunction f(nu);

        Orbital phi(SPIN::Paired);
        phi.alloc(NUMBER::Real);
        mrcpp::project(prec, phi.real(), f);
        Phi.push_back(phi);
    }


    // reference values for harmonic oscillator eigenfunctions
    DoubleMatrix ref(nFuncs,nFuncs);
    for (int i = 0; i < nFuncs; i++) {
        for (int j = 0; j < nFuncs; j++) {
            ref(i,j) = 0.0;
            if (i == j+1) ref(i,j) = std::sqrt(i/2.0);
            if (i == j-1) ref(i,j) = std::sqrt(j/2.0);
        }
    }
    PositionOperator r;
    r.setup(prec);
    SECTION("apply") {
        Orbital phi_x = r[0](Phi[0]);
        ComplexDouble X_10 = orbital::dot(Phi[1], phi_x);
        ComplexDouble X_20 = orbital::dot(Phi[2], phi_x);
        REQUIRE( X_10.real() == Approx(ref(1,0)) );
        REQUIRE( X_20.real() < thrs );
        phi_x.free();
    }
    SECTION("vector apply") {
        OrbitalVector xPhi = r[0](Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j < xPhi.size(); j++) {
                ComplexDouble X_ij = orbital::dot(Phi[i], xPhi[j]);
                REQUIRE( std::abs(X_ij.real() - ref(i,j)) < thrs);
            }
        }
        free(xPhi);
    }
    SECTION("expectation value") {
        ComplexDouble X_10 = r[0](Phi[1], Phi[0]);
        REQUIRE( X_10.real() == Approx(ref(1,0)) );
    }
    SECTION("expectation matrix ") {
        ComplexMatrix X = r[0](Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j < Phi.size(); j++) {
                REQUIRE( std::abs(X(i,j).real() - ref(i,j)) < thrs);
            }
        }
    }
    r.clear();
    free(Phi);
}

} //namespace position_operator
