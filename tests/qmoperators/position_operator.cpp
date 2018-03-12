#include "catch.hpp"

#include "mrchem.h"

#include "PositionOperator.h"
#include "Orbital.h"

using namespace mrchem;
using namespace orbital;

namespace position_operator {

// Harmonic oscillator functions

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
double N(int n) {
    double tmp = pow(2.0, n)*factorial(n)*sqrt(mrcpp::pi);
    return pow(tmp, -1.0/2.0);
}
double H(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return 2.0*x;
    return 2.0*x*H(n-1, x) - 2.0*(n-1.0)*H(n-2, x);
}

double f_0(double x) { return N(0)*H(0,x)*exp(-x*x/2.0); };
double f_1(double x) { return N(1)*H(1,x)*exp(-x*x/2.0); };
double f_2(double x) { return N(2)*H(2,x)*exp(-x*x/2.0); };
double f_3(double x) { return N(3)*H(3,x)*exp(-x*x/2.0); };

auto f_000 = [] (const double *r) -> double { return f_0(r[0])*f_0(r[1])*f_0(r[2]); };
auto f_100 = [] (const double *r) -> double { return f_1(r[0])*f_0(r[1])*f_0(r[2]); };
auto f_200 = [] (const double *r) -> double { return f_2(r[0])*f_0(r[1])*f_0(r[2]); };
auto f_300 = [] (const double *r) -> double { return f_3(r[0])*f_0(r[1])*f_0(r[2]); };

TEST_CASE("PositionOperator", "[position_operator]") {
    const double prec = 1.0e-3;
    const double thrs = prec*prec;

    OrbitalVector Phi;
    Phi.push_back(SPIN::Paired);
    Phi.push_back(SPIN::Paired);
    Phi.push_back(SPIN::Paired);
    Phi.push_back(SPIN::Paired);

    Phi[0].alloc(NUMBER::Real);
    Phi[1].alloc(NUMBER::Real);
    Phi[2].alloc(NUMBER::Real);
    Phi[3].alloc(NUMBER::Real);

    mrcpp::project(prec, Phi[0].real(), f_000);
    mrcpp::project(prec, Phi[1].real(), f_100);
    mrcpp::project(prec, Phi[2].real(), f_200);
    mrcpp::project(prec, Phi[3].real(), f_300);

    // reference values for harmonic oscillator eigenfunctions
    DoubleMatrix ref(4,4);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ref(i,j) = 0.0;
            if (i == j+1) ref(i,j) = sqrt(i/2.0);
            if (i == j-1) ref(i,j) = sqrt(j/2.0);
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
                REQUIRE( abs(X_ij.real() - ref(i,j)) < thrs);
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
                REQUIRE( abs(X(i,j).real() - ref(i,j)) < thrs);
            }
        }
    }
    r.clear();
    free(Phi);
}

} //namespace position_operator
