#include "catch.hpp"

#include "mrchem.h"

#include "IdentityOperator.h"
#include "Orbital.h"
#include "OrbitalVector.h"

using namespace mrchem;
using namespace orbital;

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

    SECTION("setup") {
        IdentityOperator I;
        REQUIRE( I.getApplyPrec() < 0.0 );

        I.setup(prec);
        REQUIRE( I.getApplyPrec() == Approx(prec) );

        I.clear();
        REQUIRE( I.getApplyPrec() < 0.0 );
    }

    SECTION("apply") {
        Orbital phi(SPIN::Paired);
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

    SECTION("vector apply") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        Phi[0].alloc(NUMBER::Real);
        Phi[1].alloc(NUMBER::Real);
        mrcpp::project(prec, Phi[0].real(), f);
        mrcpp::project(prec, Phi[1].real(), g);
        normalize(Phi);

        IdentityOperator I;
        I.setup(prec);
        SECTION("O(Phi)") {
            OrbitalVector IPhi = I(Phi);
            REQUIRE( IPhi[0].real().integrate() == Approx(Phi[0].real().integrate()) );
            REQUIRE( IPhi[1].real().integrate() == Approx(Phi[1].real().integrate()) );
            free(IPhi);
        }
        SECTION("O.dagger(Phi)") {
            OrbitalVector IPhi = I.dagger(Phi);
            REQUIRE( IPhi[0].real().integrate() == Approx(Phi[0].real().integrate()) );
            REQUIRE( IPhi[1].real().integrate() == Approx(Phi[1].real().integrate()) );
            free(IPhi);
        }
        SECTION("trace") {
            double nEl = get_electron_number(Phi);
            ComplexDouble trace = I.trace(Phi);
            REQUIRE( trace.real() == Approx(nEl) );
            REQUIRE( abs(trace.imag()) < mrcpp::MachineZero );
        }
        I.clear();
        free(Phi);
    }

    SECTION("expectation value") {
        Orbital phi(SPIN::Paired);
        phi.alloc();
        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        IdentityOperator I;
        I.setup(prec);
        SECTION("<phi|O|phi>") {
            ComplexDouble S = I(phi, phi);
            REQUIRE( S.real() == Approx(phi.squaredNorm()) );
            REQUIRE( S.imag() < mrcpp::MachineZero );
        }
        SECTION("<phi|O.dagger()|phi>") {
            ComplexDouble S = I.dagger(phi, phi);
            REQUIRE( S.real() == Approx(phi.squaredNorm()) );
            REQUIRE( S.imag() < mrcpp::MachineZero );
        }
        I.clear();
        phi.free();
    }

    SECTION("expectation matrix") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        Phi[0].alloc(NUMBER::Imag);
        Phi[1].alloc(NUMBER::Imag);
        mrcpp::project(prec, Phi[0].imag(), f);
        mrcpp::project(prec, Phi[1].imag(), g);

        IdentityOperator I;
        I.setup(prec);
        SECTION("<phi_i|O|phi_j>") {
            ComplexMatrix S = I(Phi, Phi);
            REQUIRE( abs(S(0,0)) == Approx(Phi[0].squaredNorm()) );
            REQUIRE( abs(S(1,1)) == Approx(Phi[1].squaredNorm()) );
            REQUIRE( S(0,1).real() == Approx( S(1,0).real()) );
            REQUIRE( S(0,1).imag() == Approx(-S(1,0).imag()) );
        }
        I.clear();
        free(Phi);
    }
}
