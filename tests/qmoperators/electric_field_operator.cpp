#include "catch.hpp"

#include "mrchem.h"

#include "HydrogenFunction.h"
#include "ElectricFieldOperator.h"
#include "Molecule.h"
#include "Orbital.h"
#include "orbital_utils.h"

using namespace mrchem;
using namespace orbital;

namespace electric_field_operator {

using QuantumNumbers = std::tuple<int, int, int>;
    
TEST_CASE("ElectricFieldOperator", "[electric_field_operator]") {
    const double prec = 1.0e-4;
    const double thrs = prec * prec * 10.0;
    OrbitalVector Phi;
    
    std::vector<QuantumNumbers> qn;
    qn.push_back( QuantumNumbers (1, 0, 0));
    qn.push_back( QuantumNumbers (2, 0, 0));
    qn.push_back( QuantumNumbers (2, 1, 0));
    qn.push_back( QuantumNumbers (2, 1, 1));
    qn.push_back( QuantumNumbers (2, 1, 2));
    int nFuncs = qn.size();

    Eigen::MatrixXd ref(nFuncs,nFuncs);

    //reference values for the electric field energy operator
    ref << 0.9,        0.0,  0.74493554, 0.74493554, 0.74493554,
           0.0,        0.9,  3.0,        3.0,        3.0,
           0.74493554, 3.0,  0.9,        0.0,        0.0,
           0.74493554, 3.0,  0.0,        0.9,        0.0,
           0.74493554, 3.0,  0.0,        0.0,        0.9;

    // setting up the field
    Eigen::Vector3d field;
    field << 1.0, 1.0, 1.0;
    ElectricFieldOperator EF(field);
    EF.setup(prec);

    //origin 
    double *o = new double[3];
    o[0] = 0.4;
    o[1] = 0.3;
    o[2] = 0.2;

    //setting up the orbitals
    for (int i = 0; i < nFuncs; i++) {
        HydrogenFunction f(std::get<0>(qn[i]), std::get<1>(qn[i]), std::get<2>(qn[i]), 1.0, o);
        Orbital phi(SPIN::Paired);
        phi.alloc(NUMBER::Real);
        mrcpp::project(prec, phi.real(), f);
        Phi.push_back(phi);
    }
    
    SECTION("apply") {
        Orbital phi_x = EF(Phi[0]);
        ComplexDouble X_00 = orbital::dot(Phi[0], phi_x);
        ComplexDouble X_10 = orbital::dot(Phi[1], phi_x);
        ComplexDouble X_20 = orbital::dot(Phi[2], phi_x);
        ComplexDouble X_30 = orbital::dot(Phi[3], phi_x);
        ComplexDouble X_40 = orbital::dot(Phi[4], phi_x);
        REQUIRE( X_00.real() == Approx(ref(0,0)).margin(thrs));
        REQUIRE( X_10.real() == Approx(ref(0,1)).margin(thrs));
        REQUIRE( X_20.real() == Approx(ref(0,2)).margin(thrs));
        REQUIRE( X_30.real() == Approx(ref(0,3)).margin(thrs));
        REQUIRE( X_40.real() == Approx(ref(0,4)).margin(thrs));
        phi_x.free();
    }

    SECTION("vector apply") {
        OrbitalVector xPhi = EF(Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j < xPhi.size(); j++) {
                ComplexDouble X_ij = orbital::dot(Phi[i], xPhi[j]);
                REQUIRE( std::abs(X_ij.real()) == Approx(ref(i,j)).margin(thrs));
            }
        }
        free(xPhi);
    }
    SECTION("expectation value") {
        ComplexDouble X_00 = EF(Phi[0], Phi[0]);
        REQUIRE( X_00.real() == Approx(ref(0,0)) );
    }
    SECTION("operator matrix elements") {
        ComplexMatrix X = EF(Phi, Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j < Phi.size(); j++) {
                REQUIRE( std::abs(X(i,j).real()) == Approx(ref(i,j)).margin(thrs));
            }
        }
    }
    EF.clear();
    free(Phi);
}

TEST_CASE("ElectricFieldEnergy", "[electric_field_energy]") {
    const double prec = 1.0e-4;
    const double thrs = prec * prec * 10.0;
    OrbitalVector Phi;

    // Setting up the field
    Eigen::Vector3d field;
    field << 1.0, 1.0, 1.0;
    ElectricFieldOperator EF(field);
    EF.setup(prec);

    // Setting up the molecule
    std::vector<std::string> nuc_coor;
    nuc_coor.push_back("H    1.0,    0.0,   0.0");
    nuc_coor.push_back("Li  -1.0,    0.0,   0.0");
    Molecule LiH(nuc_coor);

    double *o = new double[3];
    
    // Setting up the 1s orbital on H
    o[0] = 1.0;
    o[1] = 0.0;
    o[2] = 0.0;
    HydrogenFunction sh(1, 0, 0, 1.0, o);
    Orbital phi_h(SPIN::Paired);
    phi_h.alloc(NUMBER::Real);
    mrcpp::project(prec, phi_h.real(), sh);

    // Setting up the 1s orbital on Li
    o[0] = -1.0;
    o[1] =  0.0;
    o[2] =  0.0;
    HydrogenFunction sli(1, 0, 0, 3.0, o);
    Orbital phi_li(SPIN::Paired);
    phi_li.alloc(NUMBER::Real);
    mrcpp::project(prec, phi_li.real(), sli);

    SECTION("energy in the external field") {
        OrbitalVector Phi_LiH;
        Phi_LiH.push_back(phi_h);
        Phi_LiH.push_back(phi_li);
        double E_ext = EF.trace(Phi_LiH).real();
        double E_nex = EF.trace(LiH.getNuclei()).real();
        REQUIRE(E_ext < thrs);
        REQUIRE(E_nex == Approx(2.0));
        free(Phi_LiH);
    }
    EF.clear();
}

} //namespace electric_field_operator
