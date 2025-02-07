/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "catch2/catch_all.hpp"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>

#include "mrchem.h"

#include "analyticfunctions/HarmonicOscillatorFunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/IdentityOperator.h"
#include "qmoperators/one_electron/MomentumOperator.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include "qmoperators/one_electron/PositionOperator.h"
#include "qmoperators/one_electron/SpinOperator.h"
#include "tensor/RankTwoOperator.h"
#include "tensor/tensor_utils.h"

using namespace mrchem;

namespace operator_composition {

TEST_CASE("Operator composition", "[operator_composition]") {
    const double prec = 1.0e-3;
    const double thrs = prec * prec;

    OrbitalVector Phi;
    Phi.push_back(Orbital{SPIN::Alpha});
    Phi.push_back(Orbital{SPIN::Alpha});
    Phi.push_back(Orbital{SPIN::Alpha});
    Phi.distribute();
    for (int i = 0; i < 3; i++) {
        int nu[3] = {i, 0, 0};
        HarmonicOscillatorFunction f(nu);
        if (mrcpp::mpi::my_orb(Phi[i])) mrcpp::cplxfunc::project(Phi[i], f, NUMBER::Real, prec);
    }

    SECTION("identity operator") {
        IdentityOperator I;
        SECTION("product I*I") {
            RankZeroOperator II = I * I;
            REQUIRE(II.size() == 1);
            REQUIRE(II.size(0) == 2);

            II.setup(prec);
            const ComplexMatrix ref = I(Phi, Phi);
            const ComplexMatrix val = II(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref(0, 0).real()).epsilon(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref(0, 0).imag()).margin(thrs));
            II.clear();
        }
        SECTION("composition I(I)") {
            RankZeroOperator II = I(I);
            REQUIRE(II.size() == 1);
            REQUIRE(II.size(0) == 1);

            II.setup(prec);
            const ComplexMatrix ref = I(Phi, Phi);
            const ComplexMatrix val = II(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref(0, 0).real()).epsilon(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref(0, 0).imag()).margin(thrs));
            II.clear();
        }
        SECTION("product I*V") {
            PositionOperator V;
            RankZeroOperator IV = I * V[0];
            REQUIRE(IV.size() == 1);
            REQUIRE(IV.size(0) == 2);

            IV.setup(prec);
            const ComplexMatrix ref = V[0](Phi, Phi);
            const ComplexMatrix val = IV(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref(1, 0).real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref(1, 0).imag()).margin(thrs));
            IV.clear();
        }
        SECTION("composition I(V)") {
            PositionOperator V;
            RankZeroOperator IV = I(V[0]);
            REQUIRE(IV.size() == 1);
            REQUIRE(IV.size(0) == 1);

            IV.setup(prec);
            const ComplexMatrix ref = V[0](Phi, Phi);
            const ComplexMatrix val = IV(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref(1, 0).real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref(1, 0).imag()).margin(thrs));
            IV.clear();
        }
        SECTION("product I*D") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator ID = I * D[0];
            REQUIRE(ID.size() == 1);
            REQUIRE(ID.size(0) == 2);

            ID.setup(prec);
            const ComplexMatrix ref = D[0](Phi, Phi);
            const ComplexMatrix val = ID(Phi, Phi);
            REQUIRE(val(0, 1).real() == Catch::Approx(ref(0, 1).real()).margin(thrs));
            REQUIRE(val(0, 1).imag() == Catch::Approx(ref(0, 1).imag()).epsilon(thrs));
            ID.clear();
        }
        SECTION("composition I(D)") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator ID = I(D[0]);
            REQUIRE(ID.size() == 1);
            REQUIRE(ID.size(0) == 1);

            ID.setup(prec);
            const ComplexMatrix ref = D[0](Phi, Phi);
            const ComplexMatrix val = ID(Phi, Phi);
            REQUIRE(val(0, 1).real() == Catch::Approx(ref(0, 1).real()).margin(thrs));
            REQUIRE(val(0, 1).imag() == Catch::Approx(ref(0, 1).imag()).epsilon(thrs));
            ID.clear();
        }
        SECTION("product I*S") {
            SpinOperator S;
            RankZeroOperator IS = I * S[1];
            REQUIRE(IS.size() == 1);
            REQUIRE(IS.size(0) == 2);

            IS.setup(prec);
            const ComplexMatrix ref = S[1](Phi, Phi);
            const ComplexMatrix val = IS(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref(0, 0).real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref(0, 0).imag()).epsilon(thrs));
            IS.clear();
        }
        SECTION("composition I(S)") {
            SpinOperator S;
            RankZeroOperator IS = I(S[1]);
            REQUIRE(IS.size() == 1);
            REQUIRE(IS.size(0) == 1);

            IS.setup(prec);
            const ComplexMatrix ref = S[1](Phi, Phi);
            const ComplexMatrix val = IS(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref(0, 0).real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref(0, 0).imag()).epsilon(thrs));
            IS.clear();
        }
    }
    SECTION("potential operator") {
        PositionOperator V;
        SECTION("product V*I") {
            IdentityOperator I;
            RankZeroOperator VI = V[0] * I;
            REQUIRE(VI.size() == 1);
            REQUIRE(VI.size(0) == 2);

            VI.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexMatrix val = VI(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            VI.clear();
        }
        SECTION("composition V(I)") {
            IdentityOperator I;
            RankZeroOperator VI = V[0](I);
            REQUIRE(VI.size() == 1);
            REQUIRE(VI.size(0) == 1);

            VI.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexMatrix val = VI(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            VI.clear();
        }
        SECTION("product V*V") {
            RankZeroOperator VV = V[0] * V[0];
            REQUIRE(VV.size() == 1);
            REQUIRE(VV.size(0) == 2);

            VV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexMatrix val = VV(Phi, Phi);
            REQUIRE(val(2, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(2, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            VV.clear();
        }
        SECTION("composition V(V)") {
            RankZeroOperator VV = V[0](V[0]);
            REQUIRE(VV.size() == 1);
            REQUIRE(VV.size(0) == 1);

            VV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexMatrix val = VV(Phi, Phi);
            REQUIRE(val(2, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(2, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            VV.clear();
        }
        SECTION("product V*D") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator VD = V[0] * D[0];
            REQUIRE(VD.size() == 1);
            REQUIRE(VD.size(0) == 2);

            VD.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexMatrix val = VD(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            VD.clear();
        }
        SECTION("composition V(D)") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator VD = V[0](D[0]);
            REQUIRE(VD.size() == 1);
            REQUIRE(VD.size(0) == 2);

            VD.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexMatrix val = VD(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            VD.clear();
        }
        SECTION("product V*S") {
            SpinOperator S;
            RankZeroOperator VS = V[0] * S[2];
            REQUIRE(VS.size() == 1);
            REQUIRE(VS.size(0) == 2);

            VS.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexMatrix val = VS(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            VS.clear();
        }
        SECTION("composition V(S)") {
            SpinOperator S;
            RankZeroOperator VS = V[0](S[2]);
            REQUIRE(VS.size() == 1);
            REQUIRE(VS.size(0) == 2);

            VS.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexMatrix val = VS(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            VS.clear();
        }
    }
    SECTION("derivative operator") {
        MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
        SECTION("product D*I") {
            IdentityOperator I;
            RankZeroOperator DI = D[0] * I;
            REQUIRE(DI.size() == 1);
            REQUIRE(DI.size(0) == 2);

            DI.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 2.0};
            const ComplexMatrix val = DI(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            DI.clear();
        }
        SECTION("composition D(I)") {
            IdentityOperator I;
            RankZeroOperator DI = D[0](I);
            REQUIRE(DI.size() == 1);
            REQUIRE(DI.size(0) == 1);

            DI.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 2.0};
            const ComplexMatrix val = DI(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            DI.clear();
        }
        SECTION("product D*V") {
            PositionOperator V;
            RankZeroOperator DV = D[0] * V[0];
            REQUIRE(DV.size() == 1);
            REQUIRE(DV.size(0) == 2);

            DV.setup(prec);
            const ComplexDouble ref = {0.0, -0.5};
            const ComplexMatrix val = DV(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            DV.clear();
        }
        SECTION("composition D(V)") {
            PositionOperator V;
            RankZeroOperator DV = D[0](V[0]);
            REQUIRE(DV.size() == 1);
            REQUIRE(DV.size(0) == 1);

            DV.setup(prec);
            const ComplexDouble ref = {0.0, -1.0};
            const ComplexMatrix val = DV(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            DV.clear();
        }
        SECTION("product D*D") {
            RankZeroOperator DD = D[0] * D[0];
            REQUIRE(DD.size() == 1);
            REQUIRE(DD.size(0) == 2);

            DD.setup(prec);
            const ComplexDouble ref = {0.5, 0.0};
            const ComplexMatrix val = DD(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            DD.clear();
        }
        SECTION("composition D(D)") {
            RankZeroOperator DD = D[0](D[0]);
            REQUIRE(DD.size() == 1);
            REQUIRE(DD.size(0) == 2);

            DD.setup(prec);
            const ComplexDouble ref = {0.5, 0.0};
            const ComplexMatrix val = DD(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            DD.clear();
        }
        SECTION("product D*S") {
            SpinOperator S;
            RankZeroOperator DS = D[0] * S[2];
            REQUIRE(DS.size() == 1);
            REQUIRE(DS.size(0) == 2);

            DS.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 4.0};
            const ComplexMatrix val = DS(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            DS.clear();
        }
        SECTION("composition D(S)") {
            SpinOperator S;
            RankZeroOperator DS = D[0](S[2]);
            REQUIRE(DS.size() == 1);
            REQUIRE(DS.size(0) == 2);

            DS.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 4.0};
            const ComplexMatrix val = DS(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            DS.clear();
        }
    }
    SECTION("spin operator") {
        SpinOperator S;
        SECTION("product S*I") {
            IdentityOperator I;
            RankZeroOperator SI = S[1] * I;
            REQUIRE(SI.size() == 1);
            REQUIRE(SI.size(0) == 2);

            SI.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexMatrix val = SI(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            SI.clear();
        }
        SECTION("composition S(I)") {
            IdentityOperator I;
            RankZeroOperator SI = S[1](I);
            REQUIRE(SI.size() == 1);
            REQUIRE(SI.size(0) == 1);

            SI.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexMatrix val = SI(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            SI.clear();
        }
        SECTION("product S*V") {
            PositionOperator V;
            RankZeroOperator SV = S[2] * V[0];
            REQUIRE(SV.size() == 1);
            REQUIRE(SV.size(0) == 2);

            SV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexMatrix val = SV(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            SV.clear();
        }
        SECTION("product S(V)") {
            PositionOperator V;
            RankZeroOperator SV = S[2](V[0]);
            REQUIRE(SV.size() == 1);
            REQUIRE(SV.size(0) == 2);

            SV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexMatrix val = SV(Phi, Phi);
            REQUIRE(val(1, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(1, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            SV.clear();
        }
        SECTION("product S*D") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator SD = S[0] * D[0];
            REQUIRE(SD.size() == 1);
            REQUIRE(SD.size(0) == 2);

            SD.setup(prec);
            const ComplexDouble ref = {0.0, -std::sqrt(2.0) / 4.0};
            const ComplexMatrix val = SD(Phi, Phi);
            REQUIRE(val(0, 1).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 1).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            SD.clear();
        }
        SECTION("composition S(D)") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator SD = S[0](D[0]);
            REQUIRE(SD.size() == 1);
            REQUIRE(SD.size(0) == 2);

            SD.setup(prec);
            const ComplexDouble ref = {0.0, -std::sqrt(2.0) / 4.0};
            const ComplexMatrix val = SD(Phi, Phi);
            REQUIRE(val(0, 1).real() == Catch::Approx(ref.real()).margin(thrs));
            REQUIRE(val(0, 1).imag() == Catch::Approx(ref.imag()).epsilon(thrs));
            SD.clear();
        }
        SECTION("product S*S") {
            SpinOperator S;
            RankZeroOperator SS = S[1] * S[1];
            REQUIRE(SS.size() == 1);
            REQUIRE(SS.size(0) == 2);

            SS.setup(prec);
            const ComplexDouble ref = {1.0 / 4.0, 0.0};
            const ComplexMatrix val = SS(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            SS.clear();
        }
        SECTION("composition S(S)") {
            SpinOperator S;
            RankZeroOperator SS = S[1](S[1]);
            REQUIRE(SS.size() == 1);
            REQUIRE(SS.size(0) == 2);

            SS.setup(prec);
            const ComplexDouble ref = {1.0 / 4.0, 0.0};
            const ComplexMatrix val = SS(Phi, Phi);
            REQUIRE(val(0, 0).real() == Catch::Approx(ref.real()).epsilon(thrs));
            REQUIRE(val(0, 0).imag() == Catch::Approx(ref.imag()).margin(thrs));
            SS.clear();
        }
    }
    SECTION("vector operators") {
        NablaOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
        PositionOperator V;
        SECTION("gradient Del(V_scalar)") {
            RankOneOperator<3> gradV = D(V[0]);
            REQUIRE(gradV[0].size() == 1);
            REQUIRE(gradV[1].size() == 1);
            REQUIRE(gradV[2].size() == 1);
            REQUIRE(gradV[0].size(0) == 1);
            REQUIRE(gradV[1].size(0) == 1);
            REQUIRE(gradV[2].size(0) == 1);

            gradV.setup(prec);
            if (mrcpp::mpi::my_orb(Phi[0])) {
                OrbitalVector dPhi_0 = gradV(Phi[0]);
                REQUIRE(dPhi_0.size() == 3);
                REQUIRE(dPhi_0[0].norm() == Catch::Approx(1.0).epsilon(thrs));
                REQUIRE(dPhi_0[1].norm() == Catch::Approx(0.0).margin(thrs));
                REQUIRE(dPhi_0[2].norm() == Catch::Approx(0.0).margin(thrs));
            }
            gradV.clear();
        }
        SECTION("divergence Del . V_vector") {
            RankZeroOperator divV = tensor::dot(D, V);
            REQUIRE(divV.size() == 3);
            REQUIRE(divV.size(0) == 1);
            REQUIRE(divV.size(1) == 1);
            REQUIRE(divV.size(2) == 1);

            divV.setup(prec);
            OrbitalVector Psi = divV(Phi);
            const ComplexVector ref = orbital::get_integrals(Phi);
            const ComplexVector val = orbital::get_integrals(Psi);
            REQUIRE(val(0).real() == Catch::Approx(3.0 * ref(0).real()).epsilon(thrs));
            REQUIRE(val(0).imag() == Catch::Approx(3.0 * ref(0).imag()).margin(thrs));
            divV.clear();
        }
        SECTION("curl Del x V_vector") {
            RankOneOperator<3> curlV = tensor::cross(D, V);
            curlV.setup(prec);
            for (int i = 0; i < 3; i++) {
                RankZeroOperator curlV_i = curlV[i];
                REQUIRE(curlV[i].size() == 2);
                REQUIRE(curlV[i].size(0) == 1);
                REQUIRE(curlV[i].size(1) == 1);
                OrbitalVector Psi = curlV[i](Phi);
                const DoubleVector val = orbital::get_norms(Psi);
                for (int n = 0; n < Psi.size(); n++) REQUIRE(val(n) == Catch::Approx(0.0).margin(thrs));
            }
            curlV.clear();
        }
        SECTION("jacobian Del (x) V_vector") {
            RankTwoOperator<3, 3> jacV = tensor::outer(D, V);
            jacV.setup(prec);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    RankZeroOperator jacV_ij = jacV[i][j];
                    REQUIRE(jacV_ij.size() == 1);
                    REQUIRE(jacV_ij.size(0) == 1);
                    OrbitalVector Psi = jacV_ij(Phi);
                    const DoubleVector val = orbital::get_norms(Psi);
                    for (int n = 0; n < Psi.size(); n++) {
                        if (i == j) {
                            REQUIRE(val(n) == Catch::Approx(1.0).epsilon(thrs));
                        } else {
                            REQUIRE(val(n) == Catch::Approx(0.0).margin(thrs));
                        }
                    }
                }
            }
            jacV.clear();
        }
    }
}

} // namespace operator_composition
