#include "catch.hpp"

#include "mrchem.h"
#include "Orbital.h"

using namespace mrchem;
using namespace orbital;

namespace orbital_vector_tests {

auto f1 = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-1.0*R*R);
};

auto f2 = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-2.0*R*R);
};

auto f3 = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-3.0*R*R);
};


auto f4 = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-4.0*R*R);
};

auto f5 = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-5.0*R*R);
};

auto f6 = [] (const double *r) -> double {
    double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return exp(-6.0*R*R);
};


TEST_CASE("OrbitalVector", "[orbital_vector]") {
    const double prec = 1.0e-3;
    const double thrs = 1.0e-12;

    SECTION("push_back") {
        Orbital phi_a(SPIN::Alpha);
        Orbital phi_b(SPIN::Beta);

        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(phi_b);
        Phi.push_back(phi_a);
        Phi.push_back(phi_b);

        REQUIRE( Phi.size() == 4 );
        REQUIRE( get_electron_number(Phi, SPIN::Paired) == 5 );
        REQUIRE( get_electron_number(Phi, SPIN::Alpha) == 2 );
        REQUIRE( get_electron_number(Phi, SPIN::Beta) == 3 );
        REQUIRE( size_empty(Phi) == 0 );
        REQUIRE( size_paired(Phi) == 1 );
        REQUIRE( size_alpha(Phi) == 1 );
        REQUIRE( size_beta(Phi) == 2 );
        REQUIRE( get_multiplicity(Phi) == 2 );

        Phi.clear();
        REQUIRE( Phi.size() == 0 );
    }

    SECTION("alloc") {
        Orbital phi_a(SPIN::Alpha);
        Orbital phi_b(SPIN::Beta);

        phi_a.alloc(NUMBER::Real);
        phi_b.alloc(NUMBER::Real);
        phi_a.real().setZero();
        phi_b.real().setZero();

        OrbitalVector Phi;
        SECTION("clear") {
            Phi.push_back(phi_a);
            Phi.push_back(phi_b);
            Phi.push_back(phi_a);
            Phi.clear();
            phi_a.free();
            phi_b.free();
        }
        SECTION("free") {
            Phi.push_back(phi_a);
            Phi.push_back(phi_b);
            Phi.push_back(phi_a.deepCopy());
            free(Phi);
        }
        REQUIRE( Phi.size() == 0 );
    }

    SECTION("adjoin/disjoin vectors") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Alpha);
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Beta);
        Phi.push_back(SPIN::Beta);

        OrbitalVector Phi_p = disjoin(Phi, SPIN::Paired);
        OrbitalVector Phi_a = disjoin(Phi, SPIN::Alpha);
        OrbitalVector Phi_b = disjoin(Phi, SPIN::Beta);

        REQUIRE( Phi.size() == 0 );
        REQUIRE( Phi_p.size() == 2 );
        REQUIRE( Phi_a.size() == 1 );
        REQUIRE( Phi_b.size() == 2 );

        Phi = adjoin(Phi, Phi_p);
        Phi = adjoin(Phi, Phi_a);
        Phi = adjoin(Phi, Phi_b);

        REQUIRE( Phi.size() == 5 );
        REQUIRE( Phi_p.size() == 0 );
        REQUIRE( Phi_a.size() == 0 );
        REQUIRE( Phi_b.size() == 0 );

        REQUIRE( Phi[0].spin() == SPIN::Paired );
        REQUIRE( Phi[1].spin() == SPIN::Paired );
        REQUIRE( Phi[2].spin() == SPIN::Alpha );
        REQUIRE( Phi[3].spin() == SPIN::Beta );
        REQUIRE( Phi[4].spin() == SPIN::Beta );
    }

    SECTION("copy vectors") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Alpha);

        Phi[0].alloc(NUMBER::Real);
        Phi[1].alloc(NUMBER::Imag);
        mrcpp::project(prec, Phi[0].real(), f1);
        mrcpp::project(prec, Phi[1].imag(), f2);
        normalize(Phi);

        SECTION("copy constructor") {
            OrbitalVector Psi(Phi);
            REQUIRE( get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired) );
            REQUIRE( get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha) );
            REQUIRE( get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta) );
            REQUIRE( Psi[0].norm() == Approx(1.0) );
            REQUIRE( Psi[1].norm() == Approx(1.0) );
            Psi.clear();
        }

        SECTION("default constructor plus assignment") {
            OrbitalVector Psi;
            Psi = Phi;
            REQUIRE( get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired) );
            REQUIRE( get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha) );
            REQUIRE( get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta) );
            REQUIRE( Psi[0].norm() == Approx(1.0) );
            REQUIRE( Psi[1].norm() == Approx(1.0) );
            Psi.clear();
        }

        SECTION("default constructor plus deep copy") {
            OrbitalVector Psi;
            Psi = deep_copy(Phi);
            REQUIRE( get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired) );
            REQUIRE( get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha) );
            REQUIRE( get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta) );
            REQUIRE( Psi[0].norm() == Approx(1.0) );
            REQUIRE( Psi[1].norm() == Approx(1.0) );
            free(Psi);
        }

        SECTION("assigment constructor") {
            OrbitalVector Psi = Phi;
            REQUIRE( get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired) );
            REQUIRE( get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha) );
            REQUIRE( get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta) );
            REQUIRE( Psi[0].norm() == Approx(1.0) );
            REQUIRE( Psi[1].norm() == Approx(1.0) );
            Psi.clear();
        }

        SECTION("parameter copy") {
            OrbitalVector Psi;
            Psi = param_copy(Phi);
            REQUIRE( get_electron_number(Psi, SPIN::Paired) == get_electron_number(Phi, SPIN::Paired) );
            REQUIRE( get_electron_number(Psi, SPIN::Alpha) == get_electron_number(Phi, SPIN::Alpha) );
            REQUIRE( get_electron_number(Psi, SPIN::Beta) == get_electron_number(Phi, SPIN::Beta) );
            REQUIRE( Psi[0].norm() < 0.0);
            REQUIRE( Psi[1].norm() < 0.0);
            Psi.clear();
        }

        free(Phi);
    }

    SECTION("normalization") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Alpha);

        REQUIRE( Phi[0].norm() == Approx(-1.0) );
        REQUIRE( Phi[1].norm() == Approx(-1.0) );

        Phi[0].alloc(NUMBER::Real);
        Phi[1].alloc(NUMBER::Imag);
        mrcpp::project(prec, Phi[0].real(), f1);
        mrcpp::project(prec, Phi[1].imag(), f2);

        REQUIRE( Phi[0].norm() > 0.0 );
        REQUIRE( Phi[1].norm() > 0.0 );

        normalize(Phi);

        REQUIRE( Phi[0].norm() == Approx(1.0) );
        REQUIRE( Phi[1].norm() == Approx(1.0) );

        free(Phi);
    }

    SECTION("orthogonalization") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Beta);
        Phi.push_back(SPIN::Alpha);
        Phi.push_back(SPIN::Alpha);
        Phi.push_back(SPIN::Beta);

        Phi[0].alloc(NUMBER::Real);
        Phi[1].alloc(NUMBER::Real);
        Phi[2].alloc(NUMBER::Real);
        Phi[3].alloc(NUMBER::Real);

        mrcpp::project(prec, Phi[0].real(), f1);
        mrcpp::project(prec, Phi[1].real(), f2);
        mrcpp::project(prec, Phi[2].real(), f3);
        mrcpp::project(prec, Phi[3].real(), f4);

        orthogonalize(Phi);

        SECTION("in place orthonormalize") {
            normalize(Phi);

            int nOrbs = Phi.size();
            for (int i = 0; i < nOrbs; i++) {
                for (int j = 0; j < nOrbs; j++) {
                    ComplexDouble S_ij = orbital::dot(Phi[i], Phi[j]);
                    if (i == j ) REQUIRE( std::abs(S_ij) == Approx(1.0) );
                    if (i != j ) REQUIRE( std::abs(S_ij) < thrs );
                }
            }
        }

        SECTION("vector orthogonalize") {
            OrbitalVector Psi;
            Psi.push_back(SPIN::Alpha);
            Psi.push_back(SPIN::Beta);
            Psi[0].alloc(NUMBER::Real);
            Psi[1].alloc(NUMBER::Real);
            mrcpp::project(prec, Psi[0].real(), f5);
            mrcpp::project(prec, Psi[1].real(), f6);

            orthogonalize(Psi, Phi);

            for (int i = 0; i < Psi.size(); i++) {
                for (int j = 0; j < Phi.size(); j++) {
                    ComplexDouble S_ij = orbital::dot(Psi[i], Phi[j]);
                    REQUIRE( std::abs(S_ij) < thrs );
                }
            }
            free(Psi);
        }
        free(Phi);
    }

    SECTION("addition") {
        OrbitalVector Phi_a;
        Phi_a.push_back(SPIN::Paired);
        Phi_a.push_back(SPIN::Paired);
        Phi_a[0].alloc(NUMBER::Real);
        Phi_a[1].alloc(NUMBER::Real);
        mrcpp::project(prec, Phi_a[0].real(), f1);
        mrcpp::project(prec, Phi_a[1].real(), f2);

        OrbitalVector Phi_b;
        Phi_b.push_back(SPIN::Paired);
        Phi_b.push_back(SPIN::Paired);
        Phi_b[0].alloc(NUMBER::Real);
        Phi_b[1].alloc(NUMBER::Real);
        mrcpp::project(prec, Phi_b[0].real(), f3);
        mrcpp::project(prec, Phi_b[1].real(), f4);

        double a0 = Phi_a[0].real().integrate();
        double a1 = Phi_a[1].real().integrate();
        double b0 = Phi_b[0].real().integrate();
        double b1 = Phi_b[1].real().integrate();

        OrbitalVector Phi_c = add(1.0, Phi_a, -1.0, Phi_b);
        REQUIRE( Phi_c[0].real().integrate() == Approx(a0 - b0) );
        REQUIRE( Phi_c[1].real().integrate() == Approx(a1 - b1) );

        free(Phi_a);
        free(Phi_b);
        free(Phi_c);
    }

    SECTION("orbital transformations") {
        OrbitalVector Phi;
        Phi.push_back(SPIN::Paired);
        Phi.push_back(SPIN::Paired);
        Phi[0].alloc();
        Phi[1].alloc();
        mrcpp::project(prec, Phi[0].real(), f1);
        mrcpp::project(prec, Phi[0].imag(), f2);
        mrcpp::project(prec, Phi[1].real(), f3);
        mrcpp::project(prec, Phi[1].imag(), f4);

        orthogonalize(Phi);
        normalize(Phi);

        double theta = 0.5;
        int nOrbs = Phi.size();
        ComplexMatrix U(nOrbs, nOrbs);
        U(0,0) =  cos(theta);
        U(0,1) = -sin(theta);
        U(1,0) =  sin(theta);
        U(1,1) =  cos(theta);

        SECTION("unitary transformation") {
            OrbitalVector Psi = linear_combination(U, Phi, prec);
            for (int i = 0; i < nOrbs; i++) {
                for (int j = 0; j < nOrbs; j++) {
                    ComplexDouble S_ij = orbital::dot(Psi[i], Psi[j]);
                    if (i == j ) REQUIRE( std::abs(S_ij) == Approx(1.0) );
                    if (i != j ) REQUIRE( std::abs(S_ij) < thrs );
                }
            }
            free(Psi);
        }

        free(Phi);
    }
}

} //namespace orbital_vector_tests
