#include "catch.hpp"

#include "MRCPP/MWOperators"

#include "mrchem.h"

#include "HydrogenFunction.h"
#include "NuclearPotential.h"
#include "Orbital.h"

using namespace mrchem;
using namespace orbital;

namespace nuclear_potential {

TEST_CASE("NuclearPotential", "[nuclear_potential]") {
    const double prec = 1.0e-3;
    const double thrs = prec*prec;

    int nShells = 2;
    OrbitalVector Phi;
    for (int n = 1; n <= nShells; n++) {
        int L = n;
        for (int l = 0; l < L; l++) {
            int M = 2*l+1;
            for (int m = 0; m < M; m++) {
                HydrogenFunction f(n, l, m);

                Orbital phi(SPIN::Paired);
                phi.alloc(NUMBER::Real);
                mrcpp::project(prec, phi.real(), f);
                Phi.push_back(phi);
            }
        }
    }

    // reference values for hydrogen eigenfunctions
    int i = 0;
    DoubleVector E_P(Phi.size());
    for (int n = 1; n <= nShells; n++) {
        int L = n;
        double E_n = 1.0/(2.0*n*n);     //E_n = Z^2/(2*n^2)
        for (int l = 0; l < L; l++) {
            int M = 2*l+1;
            for (int m = 0; m < M; m++) {
                E_P(i++) = -2.0*E_n;    //virial theorem: 2<E_K> = -<E_P>
            }
        }
    }

    Nuclei nucs;
    nucs.push_back("H", 0);
    NuclearPotential V(nucs, prec);

    V.setup(prec);
    SECTION("apply") {
        Orbital Vphi_0 = V(Phi[0]);
        ComplexDouble V_00 = orbital::dot(Phi[0], Vphi_0);
        REQUIRE( V_00.real() == Approx(E_P(0)).epsilon(prec) );
        REQUIRE( V_00.imag() < thrs );
        Vphi_0.free();
    }
    SECTION("vector apply") {
        OrbitalVector VPhi = V(Phi);
        for (int i = 0; i < Phi.size(); i++) {
            ComplexDouble V_ii = orbital::dot(Phi[i], VPhi[i]);
            REQUIRE( V_ii.real() == Approx(E_P(i)).epsilon(prec) );
            REQUIRE( V_ii.imag() < thrs );
        }
        free(VPhi);
    }
    SECTION("expectation value") {
        ComplexDouble V_00 = V(Phi[0], Phi[0]);
        REQUIRE( V_00.real() == Approx(E_P(0)).epsilon(prec) );
        REQUIRE( V_00.imag() < thrs );
    }
    SECTION("expectation matrix ") {
        ComplexMatrix v = V(Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            REQUIRE( v(i,i).real() == Approx(E_P(i)).epsilon(prec) );
            REQUIRE( v(i,i).imag() < thrs );
        }
    }
    V.clear();
    free(Phi);
}

} //namespace kinetic_operator
