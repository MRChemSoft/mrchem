#include "catch.hpp"

#include "mrchem.h"
#include "Orbital.h"
#include "HydrogenFunction.h"

using namespace mrchem;

namespace density_tests {

TEST_CASE("Density", "[density]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("calc density") {
        Density rho(*MRA);

        SECTION("single orbital") {
            HydrogenFunction h(1,0,0);
            Orbital phi(SPIN::Alpha);
            phi.alloc(NUMBER::Real);
            mrcpp::project(prec, phi.real(), h);

            density::compute(-1.0, rho, phi, DENSITY::Total);
            REQUIRE( rho.integrate() == Approx(1.0) );

            phi.free();
        }
        SECTION("orbital vector") {
            HydrogenFunction h_1(2,1,0);
            HydrogenFunction h_2(2,1,1);
            HydrogenFunction h_3(2,1,2);

            OrbitalVector Phi;
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Paired);
            Phi.push_back(SPIN::Alpha);
            Phi[0].alloc(NUMBER::Real);
            Phi[1].alloc(NUMBER::Real);
            Phi[2].alloc(NUMBER::Imag);
            mrcpp::project(prec, Phi[0].real(), h_1);
            mrcpp::project(prec, Phi[1].real(), h_2);
            mrcpp::project(prec, Phi[2].imag(), h_3);

            density::compute(prec, rho, Phi, DENSITY::Total);
            REQUIRE( rho.integrate() == Approx(4.0) );

            orbital::free(Phi);
        }
    }

    SECTION("calc spin density") {
        Density rho_t(*MRA);
        Density rho_s(*MRA);
        Density rho_a(*MRA);
        Density rho_b(*MRA);

        SECTION("single orbital") {
            HydrogenFunction h(1,0,0);
            Orbital phi(SPIN::Alpha);
            phi.alloc(NUMBER::Real);
            mrcpp::project(prec, phi.real(), h);

            density::compute(-1.0, rho_t, phi, DENSITY::Total);
            density::compute(-1.0, rho_s, phi, DENSITY::Spin);
            density::compute(-1.0, rho_a, phi, DENSITY::Alpha);
            density::compute(-1.0, rho_b, phi, DENSITY::Beta);

            REQUIRE( rho_t.integrate() == Approx(1.0) );
            REQUIRE( rho_s.integrate() == Approx(1.0) );
            REQUIRE( rho_a.integrate() == Approx(1.0) );
            REQUIRE( rho_b.integrate() == Approx(0.0) );

            phi.free();
        }
        SECTION("orbital vector") {
            HydrogenFunction s1(1,0,0);
            HydrogenFunction s2(2,0,0);
            HydrogenFunction px(2,1,0);
            HydrogenFunction py(2,1,1);
            HydrogenFunction pz(2,1,2);

            OrbitalVector Phi;
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Alpha);
            Phi.push_back(SPIN::Beta);
            Phi.push_back(SPIN::Beta);

            Phi[0].alloc(NUMBER::Real);
            Phi[1].alloc(NUMBER::Real);
            Phi[2].alloc(NUMBER::Imag);
            Phi[3].alloc(NUMBER::Imag);
            Phi[4].alloc(NUMBER::Imag);
            Phi[5].alloc(NUMBER::Real);
            Phi[6].alloc(NUMBER::Real);

            mrcpp::project(prec, Phi[0].real(), s1);
            mrcpp::project(prec, Phi[1].real(), s2);
            mrcpp::project(prec, Phi[2].imag(), px);
            mrcpp::project(prec, Phi[3].imag(), py);
            mrcpp::project(prec, Phi[4].imag(), pz);
            mrcpp::project(prec, Phi[5].real(), s1);
            mrcpp::project(prec, Phi[6].real(), s2);

            density::compute(prec, rho_t, Phi, DENSITY::Total);
            density::compute(prec, rho_s, Phi, DENSITY::Spin);
            density::compute(prec, rho_a, Phi, DENSITY::Alpha);
            density::compute(prec, rho_b, Phi, DENSITY::Beta);

            REQUIRE( rho_t.integrate() == Approx(7.0) );
            REQUIRE( rho_s.integrate() == Approx(3.0) );
            REQUIRE( rho_a.integrate() == Approx(5.0) );
            REQUIRE( rho_b.integrate() == Approx(2.0) );

            orbital::free(Phi);
        }
    }
}

} //namespace orbital_tests
