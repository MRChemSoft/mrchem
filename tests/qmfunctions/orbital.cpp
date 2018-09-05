#include "catch.hpp"

#include "mrchem.h"
#include "Orbital.h"
#include "orbital_utils.h"

using namespace mrchem;

namespace orbital_tests {

auto f = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-1.0*R*R);
};

auto g = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-2.0*R*R);
};

TEST_CASE("Orbital", "[orbital]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("alloc") {
        Orbital phi_1(SPIN::Paired);
        REQUIRE( not phi_1.hasReal() );
        REQUIRE( not phi_1.hasImag() );

        phi_1.alloc(NUMBER::Real);
        REQUIRE( phi_1.hasReal() );
        REQUIRE( not phi_1.hasImag() );

        phi_1.alloc(NUMBER::Imag);
        REQUIRE( phi_1.hasReal() );
        REQUIRE( phi_1.hasImag() );

        phi_1.free(NUMBER::Real);
        REQUIRE( not phi_1.hasReal() );
        REQUIRE( phi_1.hasImag() );

        phi_1.free(NUMBER::Imag);

        REQUIRE( not phi_1.hasReal() );
        REQUIRE( not phi_1.hasImag() );
    }

    SECTION("copy") {
        Orbital phi_1(SPIN::Paired);
        phi_1.alloc(NUMBER::Real);
        mrcpp::project(prec, phi_1.real(), f);

        SECTION("copy constructor") {
            Orbital phi_2(phi_1);
            REQUIRE( phi_2.occ() == phi_1.occ() );
            REQUIRE( phi_2.spin() == phi_1.spin() );
            REQUIRE( phi_2.norm() == phi_1.norm() );
            REQUIRE( phi_2.hasReal() == phi_1.hasReal() );
            REQUIRE( phi_2.hasImag() == phi_1.hasImag() );
            phi_2.clear();
        }

        SECTION("default constructor plus assignment") {
            Orbital phi_2;
            phi_2 = phi_1;
            REQUIRE( phi_2.occ() == phi_1.occ() );
            REQUIRE( phi_2.spin() == phi_1.spin() );
            REQUIRE( phi_2.norm() == phi_1.norm() );
            REQUIRE( phi_2.hasReal() == phi_1.hasReal() );
            REQUIRE( phi_2.hasImag() == phi_1.hasImag() );
            phi_2.clear();
        }

        SECTION("default constructor plus deep copy") {
            Orbital phi_2;
            phi_2 = phi_1.deepCopy();
            REQUIRE( phi_2.occ() == phi_1.occ() );
            REQUIRE( phi_2.spin() == phi_1.spin() );
            REQUIRE( phi_2.norm() == phi_1.norm() );
            REQUIRE( phi_2.hasReal() == phi_1.hasReal() );
            REQUIRE( phi_2.hasImag() == phi_1.hasImag() );
            phi_2.free();
        }

        SECTION("assigment constructor") {
            Orbital phi_2 = phi_1;
            REQUIRE( phi_2.occ() == phi_1.occ() );
            REQUIRE( phi_2.spin() == phi_1.spin() );
            REQUIRE( phi_2.norm() == phi_1.norm() );
            REQUIRE( phi_2.hasReal() == phi_1.hasReal() );
            REQUIRE( phi_2.hasImag() == phi_1.hasImag() );
            phi_2.clear();
        }

        SECTION("parameter copy") {
            Orbital phi_2;
            phi_2 = phi_1.paramCopy();
            REQUIRE( phi_2.occ() == phi_1.occ() );
            REQUIRE( phi_2.spin() == phi_1.spin() );
            REQUIRE( phi_2.norm() < 1.0);
            REQUIRE( not phi_2.hasReal());
            REQUIRE( not phi_2.hasImag());
            phi_2.clear();
        }

        phi_1.free();
    }

    SECTION("normalize") {
        Orbital phi(SPIN::Paired);
        phi.alloc();
        REQUIRE( phi.norm() == Approx(-1.0) );

        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);
        REQUIRE( phi.norm() > 1.0 );

        phi.normalize();
        REQUIRE( phi.norm() == Approx(1.0) );

        phi.free();
    }

    SECTION("rescale") {
        Orbital phi(SPIN::Paired);
        phi.alloc();

        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        const double ref_norm = phi.norm();
        const double f_int = phi.real().integrate();
        const double g_int = phi.imag().integrate();
        SECTION("imaginary unit") {
            ComplexDouble i(0.0, 1.0);
            phi.rescale(i);
            REQUIRE( phi.norm() == Approx(ref_norm) );
            REQUIRE( phi.real().integrate() == Approx(-g_int) );
            REQUIRE( phi.imag().integrate() == Approx( f_int) );
        }
        SECTION("unitary rotation") {
            double a = sin(0.5);
            double b = cos(0.5);
            ComplexDouble i(a, b);
            phi.rescale(i);
            REQUIRE( phi.norm() == Approx(ref_norm) );
            REQUIRE( phi.real().integrate() == Approx(a*f_int - b*g_int) );
            REQUIRE( phi.imag().integrate() == Approx(b*f_int + a*g_int) );
        }

        phi.free();
    }

    SECTION("orthogonalize") {
        Orbital phi_1(SPIN::Alpha);
        phi_1.alloc(NUMBER::Real);
        mrcpp::project(prec, phi_1.real(), f);

        SECTION("different spin") {
            Orbital phi_2(SPIN::Beta);
            phi_2.alloc(NUMBER::Imag);
            mrcpp::project(prec, phi_2.imag(), g);

            ComplexDouble S = orbital::dot(phi_1, phi_2);
            REQUIRE( std::abs(S.real()) < thrs );
            REQUIRE( std::abs(S.imag()) < thrs );

            phi_2.free();
        }

        SECTION("same spin") {
            Orbital phi_2(SPIN::Alpha);
            phi_2.alloc(NUMBER::Imag);
            mrcpp::project(prec, phi_2.imag(), g);

            ComplexDouble S1 = orbital::dot(phi_1, phi_2);
            REQUIRE( std::abs(S1.real()) < thrs );
            REQUIRE( std::abs(S1.imag()) > thrs );

            ComplexDouble S2 = orbital::dot(phi_1, phi_2.dagger());
            REQUIRE( S2.real() == Approx( S1.real()) );
            REQUIRE( S2.imag() == Approx(-S1.imag()) );

            phi_2.orthogonalize(phi_1);

            ComplexDouble S3 = orbital::dot(phi_1, phi_2);
            REQUIRE( std::abs(S3.real()) < thrs );
            REQUIRE( std::abs(S3.imag()) < thrs );

            phi_2.free();
        }
        phi_1.free();
    }

    SECTION("add") {
        ComplexDouble c(0.5, 0.5);
        Orbital phi(SPIN::Paired);
        phi.alloc();
        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        SECTION("scalar conjugate") {
            Orbital psi = orbital::add(std::conj(c), phi, c, phi);
            REQUIRE( psi.real().integrate() == Approx(phi.real().integrate()) );
            REQUIRE( psi.imag().integrate() == Approx(phi.imag().integrate()) );
            psi.free();
        }

        SECTION("orbital conjugate") {
            Orbital psi = orbital::add(c, phi, c, phi.dagger());
            REQUIRE( psi.real().integrate() == Approx(phi.real().integrate()) );
            REQUIRE( psi.imag().integrate() == Approx(phi.real().integrate()) );
            psi.free();
        }
        phi.free();
    }

    SECTION("multiply") {
        Orbital phi(SPIN::Paired);
        phi.alloc();
        mrcpp::project(prec, phi.real(), f);
        mrcpp::project(prec, phi.imag(), g);

        Orbital psi = orbital::multiply(phi.dagger(), phi);
        double f_norm = phi.real().getSquareNorm();
        double g_norm = phi.imag().getSquareNorm();
        REQUIRE( psi.real().integrate() == Approx(f_norm + g_norm) );
        REQUIRE( psi.imag().integrate() < thrs );

        phi.free();
        psi.free();
    }
}

} //namespace orbital_tests
