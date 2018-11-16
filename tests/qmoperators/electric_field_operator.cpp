/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "catch.hpp"

#include "mrchem.h"
#include "parallel.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using namespace mrchem;
using namespace orbital;

namespace electric_field_operator {

using QuantumNumbers = std::tuple<int, int, int>;

TEST_CASE("ElectricFieldOperator", "[electric_field_operator]") {
    const double prec = 1.0e-4;
    const double thrs = prec * prec * 10.0;
    OrbitalVector Phi;

    std::vector<QuantumNumbers> qn;
    qn.push_back(QuantumNumbers(1, 0, 0));
    qn.push_back(QuantumNumbers(2, 0, 0));
    qn.push_back(QuantumNumbers(2, 1, 0));
    qn.push_back(QuantumNumbers(2, 1, 1));
    qn.push_back(QuantumNumbers(2, 1, 2));
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
    ElectricFieldOperator EF(field, {0.0, 0.0, 0.0});
    EF.setup(prec);

    //origin
    mrcpp::Coord<3> o{0.4, 0.3, 0.2};

    //setting up the orbitals
    for (int i = 0; i < nFuncs; i++) Phi.push_back(SPIN::Paired);
    mpi::distribute(Phi);

    for (int i = 0; i < nFuncs; i++) {
        HydrogenFunction f(std::get<0>(qn[i]), std::get<1>(qn[i]), std::get<2>(qn[i]), 1.0, o);
        if (mpi::my_orb(Phi[i])) Phi[i].alloc(NUMBER::Real);
        if (mpi::my_orb(Phi[i])) mrcpp::project(prec, Phi[i].real(), f);
    }

    SECTION("apply") {
        //update ref based on MPI
        for (int i = 0; i < nFuncs; i++) {
            for (int j = 0; j < nFuncs; j++) {
                if (not mpi::my_orb(Phi[i]) or not mpi::my_orb(Phi[j])) ref(i,j) = 0.0;
            }
        }

        Orbital phi_x = EF(Phi[0]);
        ComplexDouble X_00 = orbital::dot(Phi[0], phi_x);
        ComplexDouble X_10 = orbital::dot(Phi[1], phi_x);
        ComplexDouble X_20 = orbital::dot(Phi[2], phi_x);
        ComplexDouble X_30 = orbital::dot(Phi[3], phi_x);
        ComplexDouble X_40 = orbital::dot(Phi[4], phi_x);
        REQUIRE(X_00.real() == Approx(ref(0,0)).margin(thrs));
        REQUIRE(X_10.real() == Approx(ref(0,1)).margin(thrs));
        REQUIRE(X_20.real() == Approx(ref(0,2)).margin(thrs));
        REQUIRE(X_30.real() == Approx(ref(0,3)).margin(thrs));
        REQUIRE(X_40.real() == Approx(ref(0,4)).margin(thrs));
        phi_x.free();
    }

    SECTION("vector apply") {
        //update ref based on MPI
        for (int i = 0; i < nFuncs; i++) {
            for (int j = 0; j < nFuncs; j++) {
                if (not mpi::my_orb(Phi[i]) or not mpi::my_orb(Phi[j])) ref(i,j) = 0.0;
            }
        }

        OrbitalVector xPhi = EF(Phi);
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j < xPhi.size(); j++) {
                ComplexDouble X_ij = orbital::dot(Phi[i], xPhi[j]);
                REQUIRE(std::abs(X_ij.real()) == Approx(ref(i,j)).margin(thrs));
            }
        }
        free(xPhi);
    }
    SECTION("expectation value") {
        ComplexDouble X_00 = EF(Phi[0], Phi[0]);
        if (mpi::my_orb(Phi[0])) {
            REQUIRE(X_00.real() == Approx(ref(0,0)));
        } else {
            REQUIRE(X_00.real() == Approx(0.0).margin(thrs));
        }
    }
    SECTION("operator matrix elements") {
        ComplexMatrix X = EF(Phi, Phi);
        for (int i = 0; i < X.rows(); i++) {
            for (int j = 0; j < X.cols(); j++) {
                REQUIRE(std::abs(X(i,j).real()) == Approx(ref(i,j)).margin(thrs));
            }
        }
    }
    EF.clear();
    free(Phi);
}

TEST_CASE("ElectricFieldEnergy", "[electric_field_energy]") {
    const double prec = 1.0e-4;
    const double thrs = prec * prec * 10.0;

    // Setting up the field
    Eigen::Vector3d field;
    field << 1.0, 1.0, 1.0;
    ElectricFieldOperator EF(field, {0.0, 0.0, 0.0});
    EF.setup(prec);

    // Setting up the molecule
    std::vector<std::string> nuc_coor;
    nuc_coor.push_back("H    1.0,    0.0,   0.0");
    nuc_coor.push_back("Li  -1.0,    0.0,   0.0");
    Molecule LiH(nuc_coor);

    double *o = new double[3];

    OrbitalVector Phi;
    Phi.push_back(SPIN::Paired);
    Phi.push_back(SPIN::Paired);
    mpi::distribute(Phi);

    // Setting up the 1s orbital on H
    HydrogenFunction sh(1, 0, 0, 1.0, {1.0, 0.0, 0.0});
    Orbital &phi_h = Phi[0];
    if (mpi::my_orb(phi_h)) phi_h.alloc(NUMBER::Real);
    if (mpi::my_orb(phi_h)) mrcpp::project(prec, phi_h.real(), sh);

    // Setting up the 1s orbital on Li
    HydrogenFunction sli(1, 0, 0, 3.0, {-1.0, 0.0, 0.0});
    Orbital &phi_li = Phi[1];
    if (mpi::my_orb(phi_li)) phi_li.alloc(NUMBER::Real);
    if (mpi::my_orb(phi_li)) mrcpp::project(prec, phi_li.real(), sli);

    SECTION("energy in the external field") {
        double E_ext = EF.trace(Phi).real();
        double E_nex = EF.trace(LiH.getNuclei()).real();
        REQUIRE(E_ext < thrs);
        REQUIRE(E_nex == Approx(2.0));
        free(Phi);
    }
    EF.clear();
}

} //namespace electric_field_operator
