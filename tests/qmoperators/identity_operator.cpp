#include "catch.hpp"

#include "mrchem.h"
#include "mrtest.h"

#include "IdentityOperator.h"
#include "Orbital.h"

extern MultiResolutionAnalysis<3> *mrchem::MRA;  //< Default MRA

using namespace mrchem;

auto f = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-1.0*R*R);
};

auto g = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-2.0*R*R);
};

TEST_CASE("IdentityOperator", "[identity_operator]") {
    const double prec = 1.0e-3;
    mrchem::MRA = mrtest::initialize_mra();

    SECTION("setup") {
        IdentityOperator I;
        REQUIRE( I.getApplyPrec() < 0.0 );

        I.setup(prec);
        REQUIRE( I.getApplyPrec() == Approx(prec) );

        I.clear();
        REQUIRE( I.getApplyPrec() < 0.0 );
    }

    SECTION("apply") {
        Orbital phi(2, SPIN::Paired);
        phi.alloc();
        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        IdentityOperator I;
        I.setup(prec);
        SECTION("O(phi)") {
            Orbital Iphi = I(phi);
            REQUIRE( Iphi.real().integrate() == Approx(phi.real().integrate()) );
            REQUIRE( Iphi.imag().integrate() == Approx(phi.imag().integrate()) );
            Iphi.free();
        }
        SECTION("O.dagger(phi)") {
            Orbital Iphi = I.dagger(phi);
            REQUIRE( Iphi.real().integrate() == Approx(phi.real().integrate()) );
            REQUIRE( Iphi.imag().integrate() == Approx(phi.imag().integrate()) );
            Iphi.free();
        }
        I.clear();
        phi.free();
    }

    SECTION("expectation value") {
        Orbital phi(2, SPIN::Paired);
        phi.alloc();
        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        IdentityOperator I;
        I.setup(prec);
        SECTION("<phi|O|phi>") {
            ComplexDouble S = I(phi, phi);
            REQUIRE( S.real() == Approx(phi.squaredNorm()) );
            REQUIRE( S.imag() == Approx(0.0) );
        }
        SECTION("<phi|O.dagger()|phi>") {
            ComplexDouble S = I.dagger(phi, phi);
            REQUIRE( S.real() == Approx(phi.squaredNorm()) );
            REQUIRE( S.imag() == Approx(0.0) );
        }
        I.clear();
        phi.free();
    }

    delete mrchem::MRA;
}
