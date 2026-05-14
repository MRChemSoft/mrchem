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

#include "FockBuilder.h"

#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include "MRCPP/utils/spinor_utils.h"

#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "ReactionOperator.h"
#include "XCOperator.h"
#include "analyticfunctions/NuclearFunction.h"
#include "chemistry/chemistry_utils.h"
#include "properties/SCFEnergy.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/IdentityOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/one_electron/ZoraOperator.h"
#include "qmoperators/qmoperator_utils.h"
#include "utils/math_utils.h"

#include "qmoperators/one_electron/AZoraPotential.h"

#include <filesystem>

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief build the Fock operator once all contributions are in place
 *
 */
void FockBuilder::build(double exx) {
    std::cout << "FockBuilder::build -- Building Potential operator" << std::endl;
    this->exact_exchange = exx;

    this->V = RankZeroOperator();
    if (this->nuc != nullptr) this->V += (*this->nuc); //TODO, peut-être tenter de retirer chacune de ces contributions pour voir laquelle est mal initialisée
    if (this->coul != nullptr) this->V += (*this->coul); //coulomb pose un problème
    if (this->ex != nullptr) this->V -= this->exact_exchange * (*this->ex);
    if (this->xc != nullptr) this->V += (*this->xc);
    if (this->ext != nullptr) this->V += (*this->ext);
    if (this->Ro != nullptr) this->V -= (*this->Ro);
    //debug, need to be removed eventually
    if (this->nuc == nullptr) std::cout << "FockBuilder::build -- Nuclear operator not present" << std::endl;
    if (this->coul == nullptr) std::cout << "FockBuilder::build -- Coulomb operator not present" << std::endl;
    if (this->ex == nullptr) std::cout << "FockBuilder::build -- Exchange operator not present" << std::endl;
    if (this->nuc != nullptr) std::cout << "FockBuilder::build -- Nuclear operator present" << std::endl;
    if (this->coul != nullptr) std::cout << "FockBuilder::build -- Coulomb operator present" << std::endl;
    if (this->ex != nullptr) std::cout << "FockBuilder::build -- Exchange operator present" << std::endl;
    // if (this->xc == nullptr) std::cout << "FockBuilder::build -- XC operator not present" << std::endl;
    // if (this->ext == nullptr) std::cout << "FockBuilder::build -- External operator not present" << std::endl;
    // if (this->Ro == nullptr) std::cout << "FockBuilder::build -- Reaction operator not present" << std::endl;
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will call the setup function of all underlying operators, and in particular
 * it will compute the internal exchange if there is an ExchangeOperator.
 */
void FockBuilder::setup(double prec) {
    Timer t_tot;

    std::cout << "FockBuilder::setup -- ex is nullptr =" << (this->ex== nullptr)<< "-- coul is nullptr =" << (this->coul== nullptr) << std::endl;

    // std::cout << "FockBuilder::setup -- Starting setup of kinetic and potential operators Potential=" << (this->V == nullptr) << std::endl;

    auto plevel = Printer::getPrintLevel();
    if (plevel == 2) {
        mrcpp::print::header(2, "Building Fock operator");
        mrcpp::print::value(2, "Precision", prec, "(rel)", 5);
        mrcpp::print::separator(2, '-');
    }
    // // //test TODOD À exécuter
    // for (auto &i : this->coul()->orbitals) {
    //     MSG_INFO("Orbital is ")
    //     // for (int j = 0; j < i.size(); j++) { std::cout << "FockBuilder::setup -- operator: "<< j << " " << i[j] << std::endl; }
    // }
    // std::cout << "FockBuilder::setup -- Starting setup of operators: " << &i << std::endl;
    // if (this->coul != nullptr) this->V += (*this->coul);
    // std::cout << "FockBuilder::setup -- Starting setup of kinetic and potential operators" << std::endl;
    this->prec = prec;
    // std::cout << "FockBuilder::setup -- Starting setup of kinetic and potential operators 2 "<< (this->mom != nullptr) << std::endl;
    if (this->mom != nullptr) this->momentum().setup(prec);
    // std::cout << "FockBuilder::setup -- Starting setup of kinetic and potential operators 3" << std::endl;
    this->potential().setup(prec);
    // std::cout << "FockBuilder::setup -- Starting setup of kinetic and potential operators 4" << std::endl;
    this->perturbation().setup(prec); //TODO: uncomment when response is implemented for multiple components

    // std::cout << "FockBuilder::setup -- Kinetic and potential operators setup done" << std::endl;

    if (isZora()) {
        // MSG_INFO("Setting up ZORA operators");
        Timer t_zora; //TODO: make this working for 2C
        double c = getLightSpeed();
        // MSG_INFO("c ok");
        mrcpp::print::header(3, "Building ZORA operators");
        mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
        mrcpp::print::value(3, "Light speed", c, "(au)", 5);
        mrcpp::print::separator(3, '-');
        auto vz = collectZoraBasePotential();
        // MSG_INFO("ZORA base potential collected");
        // chi = kappa - 1. See ZoraOperator.h for more information.
        this->chi = std::make_shared<ZoraOperator>(*vz, c, prec, false);
        // MSG_INFO("ZORA chi operator setup");
        this->chi_inv = std::make_shared<ZoraOperator>(*vz, c, prec, true);
        // MSG_INFO("ZORA chi inverse operator setup");
        this->zora_base = RankZeroOperator(vz);
        // MSG_INFO("ZORA base operator setup");
        this->chi->setup(prec);
        // MSG_INFO("ZORA chi setup ok");
        this->chi_inv->setup(prec);
        // MSG_INFO("ZORA chi inverse setup ok");
        this->zora_base.setup(prec);
        // MSG_INFO("ZORA base setup ok");
        mrcpp::print::footer(3, t_zora, 2);
    }
    if (isAZora()) {
        Timer t_zora;
        double c = getLightSpeed();
        mrcpp::print::header(3, "Building AZORA operators");
        mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
        mrcpp::print::value(3, "Light speed", c, "(au)", 5);
        mrcpp::print::separator(3, '-');
        int adap = 0;

        chiPot->project(prec);
        chiInvPot = std::make_shared<QMPotential>(adap);

        mrcpp::deep_copy(*chiInvPot, *chiPot);

        chiInvPot->real().map([](double val) { return 1.0 / (val + 1) - 1; });

        this->chi = std::make_shared<ZoraOperator>(chiPot, "kappa");
        this->chi_inv = std::make_shared<ZoraOperator>(chiInvPot, "kappa_inv");
        this->chi->setup(prec);
        this->chi_inv->setup(prec);

        mrcpp::print::footer(3, t_zora, 2);
    }

    t_tot.stop();
    if (plevel == 2) mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_tot);
}

/** @brief clear operator after application
 *
 * This will call the clear function of all underlying operators, and bring them back
 * to the state after construction. The operator can now be reused after another setup.
 */
void FockBuilder::clear() {
    if (this->mom != nullptr) this->momentum().clear();
    this->potential().clear();
    this->perturbation().clear();
    if (isZora()) {
        this->chi->clear();
        this->chi_inv->clear();
        this->zora_base.clear();
    }
    if (isAZora()) {
        chi->clear();
        chi_inv->clear();
        chiPot->free(mrchem::NUMBER::Total);
        chiInvPot->free(mrchem::NUMBER::Total);
    }
}

/** @brief rotate orbitals of two-electron operators
 *
 * @param U: unitary transformation matrix
 *
 * This function should be used in case the orbitals are rotated *after* the FockBuilder
 * has been setup. In particular the ExchangeOperator needs to rotate the precomputed
 * internal exchange potentials.
 */
void FockBuilder::rotate(const ComplexMatrix &U) {
    if (this->ex != nullptr) this->ex->rotate(U);
}

/** @brief compute the SCF energy
 *
 * @param Phi: orbitals
 * @param F: Fock matrix
 *
 * This function will compute the total energy for a given OrbitalVector and
 * the corresponding Fock matrix. Tracing the kinetic energy operator is avoided
 * by tracing the Fock matrix and subtracting all other contributions.
 */
SCFEnergy FockBuilder::trace(OrbitalVector &Phi, const Nuclei &nucs) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing molecular energy");

    double E_kin = 0.0;  // Kinetic energy
    double E_nn = 0.0;   // Nuclear repulsion
    double E_en = 0.0;   // Nuclear-electronic interaction
    double E_ee = 0.0;   // Electronic repulsion
    double E_x = 0.0;    // Exact Exchange
    double E_xc = 0.0;   // Exchange and Correlation
    double E_eext = 0.0; // External field contribution to the electronic energy
    double E_next = 0.0; // External field contribution to the nuclear energy
    double Er_nuc = 0.0; // Nuclear reaction energy
    double Er_el = 0.0;  // Electronic reaction energy
    double Er_tot = 0.0; // Total reaction energy

    // Nuclear part
    MSG_INFO("nuc");
    if (this->nuc != nullptr) E_nn = chemistry::compute_nuclear_repulsion(nucs);
    if (this->ext != nullptr) E_next = -this->ext->trace(nucs).real();

    // Reaction potential part
    if (this->Ro != nullptr) {
        Density rho_el(false);
        density::compute(this->prec, rho_el, Phi, DensityType::Total);
        rho_el.rescale(-1.0);
        std::tie(Er_el, Er_nuc) = this->Ro->getSolver()->computeEnergies(rho_el);

        Er_tot = Er_nuc + Er_el;
    }

    // Kinetic part
    MSG_INFO("kin");
    if (isZora() || isAZora()) {
        bool spinorial = (Phi[0].Ncomp() > 1); //assumes all orbitals have the same number of components
        //second term doesn't inclue Pauli matrices (i.e. spinorial is false) because (σ·p)(σ·p) = p^2
        E_kin = qmoperator::calc_kinetic_trace(momentum(), *this->chi, Phi, spinorial).real() + qmoperator::calc_kinetic_trace(momentum(), Phi);
    } else {
        E_kin = qmoperator::calc_kinetic_trace(momentum(), Phi);
    }

    // Electronic part
    MSG_INFO("ee");
    if (this->nuc != nullptr) { E_en = this->nuc->trace(Phi).real(); }

    if (this->coul != nullptr) E_ee = 0.5 * this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    // ComplexDouble tutex = this->ex->trace(Phi).real();
    // MSG_INFO("coulomb expct val="<< E_ee << " exchange expct val=" << -this->exact_exchange << " " << tutex.real() << tutex.imag());
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->ext != nullptr) E_eext = this->ext->trace(Phi).real();
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing molecular energy", t_tot);
    MSG_INFO("trace ok");

    return SCFEnergy{E_kin, E_nn, E_en, E_ee, E_x, E_xc, E_next, E_eext, Er_tot, Er_nuc, Er_el};
}

ComplexMatrix FockBuilder::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix");

    // MSG_INFO("pre kinetic mat");
    ComplexMatrix T_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    if (isZora() || isAZora()) {
        //If we have spinors, the kinetic operator is of the form (σ·p)V(σ·p), with σ being a Pauli matrix.
        //What this boolean does is enabling the application of the Pauli matrices along the x,y,z momentum operators.
        //NOTE! The second term does not change from being spinorial; (σ·p)(σ·p) = p^2 using the Dirac identity.
        bool spinorial = (bra[0].Ncomp() > 1); //assumes all orbitals have the same number of components
        T_mat = (qmoperator::calc_kinetic_matrix(momentum(), *this->chi, bra, ket, spinorial) + qmoperator::calc_kinetic_matrix(momentum(), bra, ket));
    } else {
        T_mat = qmoperator::calc_kinetic_matrix(momentum(), bra, ket);
    }

    //debug tests start
    for (int a = 0; a < T_mat.rows(); a++) {
        for (int b = 0; b < T_mat.cols(); b++) {
            std::cout<< "T_mat(" << a << ", " << b << ") = " << T_mat(a, b) << "; ";
        }
        std::cout << std::endl;
    }
    MSG_INFO("post kin mat, pre pot mat");
    ComplexMatrix Vnuc_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    Vnuc_mat = (*getNuclearOperator())(bra, ket);
    for (int a = 0; a < Vnuc_mat.rows(); a++) {
        for (int b = 0; b < Vnuc_mat.cols(); b++) {
            std::cout<< "Vnuc_mat(" << a << ", " << b << ") = " << Vnuc_mat(a, b) << "; ";
        }
        std::cout << std::endl;
    }
    ComplexMatrix Vcoul_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    Vcoul_mat = (*getCoulombOperator())(bra, ket);
    for (int a = 0; a < Vcoul_mat.rows(); a++) {
        for (int b = 0; b < Vcoul_mat.cols(); b++) {
            std::cout<< "Vcoul_mat(" << a << ", " << b << ") = " << Vcoul_mat(a, b) << "; ";
        }
        std::cout << std::endl;
    }
    ComplexMatrix Vex_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    Vex_mat = (*getExchangeOperator())(bra, ket);
    for (int a = 0; a < Vex_mat.rows(); a++) {
        for (int b = 0; b < Vex_mat.cols(); b++) {
            std::cout<< "Vex_mat(" << a << ", " << b << ") = " << Vex_mat(a, b) << "; ";
        }
        std::cout << std::endl;
    }
    //debug tests end

    ComplexMatrix V_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    V_mat += potential()(bra, ket);

    // debug tests start
    for (int a = 0; a < V_mat.rows(); a++) {
        for (int b = 0; b < V_mat.cols(); b++) {
            std::cout<< "V_mat(" << a << ", " << b << ") = " << V_mat(a, b) << "; ";
        }
        std::cout << std::endl;
    }
    MSG_INFO("ok");
    // debug tests end

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Fock matrix", t_tot);
    return T_mat + V_mat;
}

OrbitalVector FockBuilder::buildHelmholtzArgument(double prec, OrbitalVector Phi, ComplexMatrix F_mat, ComplexMatrix L_mat) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Helmholtz argument");

    Timer t_rot;
    MSG_INFO("Pouet")
    //Debug test
    // MSG_INFO("Fock matrix ");
    for (int a = 0; a < F_mat.rows(); a++) {
        for (int b = 0; b < F_mat.cols(); b++) {
            std::cout<< "F_mat(" << a << ", " << b << ") = " << F_mat(a, b) << "; ";
        }
        std::cout << std::endl;
    }
    // MSG_INFO("Overlap? matrix ");
    for (int a = 0; a < L_mat.rows(); a++) {
        for (int b = 0; b < L_mat.cols(); b++) {
            std::cout<< "L_mat(" << a << ", " << b << ") = " << L_mat(a, b) << "; ";
        }
        std::cout << std::endl;
    }
    //end debug test
    OrbitalVector Psi = orbital::rotate(Phi, L_mat - F_mat, prec);
    // MSG_INFO("pouetronicus")
    mrcpp::print::time(2, "Rotating orbitals", t_rot);

    OrbitalVector out;
    if (isZora() || isAZora()) {
        out = buildHelmholtzArgumentZORA(Phi, Psi, F_mat.real().diagonal(), prec);
    } else {
        out = buildHelmholtzArgumentNREL(Phi, Psi);
    }
    Psi.clear();

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Helmholtz argument", t_tot);
    return out;
}

/**
 * @brief Build the Helmholtz argument for the ZORA operator. Eq. 17 in J. Chem. Theory and Comput. 2024, 20, 728-737
 */
OrbitalVector FockBuilder::buildHelmholtzArgumentZORA(OrbitalVector &Phi, OrbitalVector &Psi, DoubleVector eps, double prec) {
    // Get necessary operators
    double c = getLightSpeed();
    double two_cc = 2.0 * c * c;
    MomentumOperator &p = momentum();
    RankZeroOperator &V = potential();
    RankZeroOperator &chi = *this->chi;
    RankZeroOperator &chi_m1 = *this->chi_inv;
    RankZeroOperator operOne = 0.5 * tensor::dot(p(chi), p);
    // MSG_INFO("start");

    std::shared_ptr<RankZeroOperator> operThreePtr = nullptr;

    if (isZora()) {
        RankZeroOperator &V_zora = this->zora_base;
        // MSG_INFO("V_zora initialised");
        operThreePtr = std::make_shared<RankZeroOperator>(V_zora * chi + V_zora);
        // MSG_INFO("V_zora computed");
    } else if (isAZora()) {
        /*
        Note that V_z * kappa = 2 c^2 * (kappa - 1)
        With this trick, the expensive projection of the potential is avoided
        */
        std::shared_ptr<QMPotential> vTimesKappa = std::make_shared<QMPotential>(0);
        mrcpp::deep_copy(*vTimesKappa, *chiPot);
        vTimesKappa->real().map([two_cc](double val) { return two_cc * (val); });
        operThreePtr = std::make_shared<RankZeroOperator>(vTimesKappa);
    } else {
        MSG_ABORT("At this point, the ZORA or AZORA operator should be set. Exiting.");
    }

    RankZeroOperator operThree = *operThreePtr;

    operOne.setup(prec);
    operThree.setup(prec);
    // MSG_INFO("1 & 3 setup");

    // Compute OrbitalVectors
    Timer t_1;
    OrbitalVector termOne = operOne(Phi);
    // MSG_INFO("1 applied");
    // for (int i = 0; i < termOne.size(); i++) {
    //     //termOne will have a prefactor -1, because of i*i
    //     if (not mrcpp::mpi::my_func(termOne[i])) continue;
    //     ComplexDouble fac = termOne[i].func_ptr->data.c1[0];
    //     if(std::norm(fac-1.0)>mrcpp::MachineZero)termOne[i].rescale(fac);
    //     termOne[i].func_ptr->data.c1[0] = {1.0, 0.0};
    // }


    mrcpp::print::time(2, "Computing gradient term", t_1);

    Timer t_2;
    OrbitalVector termTwo = V(Phi);
    // MSG_INFO("2 applied");

    mrcpp::print::time(2, "Computing potential term", t_2);

    // Compute transformed orbitals scaled by diagonal Fock elements
    Timer t_3;
    OrbitalVector epsPhi = orbital::deep_copy(Phi);
    for (int i = 0; i < epsPhi.size(); i++) {
        if (not mrcpp::mpi::my_func(epsPhi[i])) continue;
        epsPhi[i].rescale(eps[i] / two_cc);
    }
    OrbitalVector termThree = operThree(epsPhi);
    // MSG_INFO("3 applied");
    mrcpp::print::time(2, "Computing rescaled potential term", t_3);

    //spin orbit coupling term, which would be identically 0 for scalar functions.
    OrbitalVector termSO(Phi.size());
    // OrbitalVector termSO = orbital::deep_copy(Phi); //test debug test start
    //set the orbitals to zero in case we don't use them 
    // for (int i = 0; i < termSO.size(); i++) {
    //     //def termSO real or complex depending on phi ? 
    //     // termSO[i].alloc(Phi[i].Ncomp(), true);
    //     for (int comp = 0; comp < termSO[i].Ncomp(); comp++){
    //         if (termSO[i].isreal()){termSO[i].CompD[comp]->setZero();}
    //         if (termSO[i].iscomplex()){termSO[i].CompC[comp]->setZero();}
    //     }
    // }
    //test debug test end
    if ((Phi[0].Ncomp() == 2) and isZora()) {
        // mrcpp::apply(prec, *dx_chi, p[0], chi);
        
        // MSG_INFO("spinorial");
        //Manually implementing the cross product appearing in the spin-orbit term
        //NOTE! The multiplication by the Pauli matrices will need to be handled later,
        //      during the application of the cross-product to the orbitals.
        // auto dchi = p(chi); 
        RankZeroOperator operSOX = p(chi)[1]*p[2] - p(chi)[2]*p[1];
        // RankZeroOperator operSOX = dchi[1]*p[2] - dchi[2]*p[1];
        operSOX.setup(prec);
        // MSG_INFO("Spinorbit =" << " " << (p(chi)[0].getOperatorExpansion()[0])->getSquareNorm()); //<< (p(chi)[1]).getSquareNorm()<< (p(chi)[2]).getSquareNorm());
        RankZeroOperator operSOY = p(chi)[2]*p[0] - p(chi)[0]*p[2];
        // RankZeroOperator operSOY = dchi[2]*p[0] - dchi[0]*p[2];
        operSOY.setup(prec);
        RankZeroOperator operSOZ = p(chi)[0]*p[1] - p(chi)[1]*p[0];
        // RankZeroOperator operSOZ = dchi[0]*p[1] - dchi[1]*p[0];
        operSOZ.setup(prec);
        //Applying the curls to temporary copies of the orbitals 
        //NOTE! Extremely inefficient!
        // OrbitalVector orbTempX = orbital::deep_copy(Phi);
        OrbitalVector orbTempX = operSOX(Phi);
        // orbTempX = operSOX(Phi);
        // OrbitalVector orbTempY = orbital::deep_copy(Phi);
        OrbitalVector orbTempY = operSOY(Phi);
        // orbTempY = operSOY(Phi);
        // OrbitalVector orbTempZ = orbital::deep_copy(Phi);
        OrbitalVector orbTempZ = operSOZ(Phi);
        // orbTempZ = operSOZ(Phi);
        // MSG_INFO("Curls applied");
        //adding the contributions together. Note that the coefficient is purely imaginary.
        for (int i = 0; i < Phi.size(); i++) {
            // termSO[i].defcomplex(); //test debug test
            // termSO[i].alloc(Phi[i].Ncomp(), true);

            Orbital orbTempX2;//test debug test
            Orbital orbTempY2;//test debug test
            Orbital orbTempZ2;//test debug test
            deep_copy(orbTempX2, orbTempX[i]);//test debug test
            deep_copy(orbTempY2, orbTempY[i]);//test debug test
            deep_copy(orbTempZ2, orbTempZ[i]);//test debug test
            // mrcpp::apply_Pauli(orbTempX[i], orbTempX[i], 1, -1.0, false);//it runs but I believe it is not computing correctly if input and output is the same
            // mrcpp::apply_Pauli(orbTempY[i], orbTempY[i], 1, -1.0, false); //runs
            // mrcpp::apply_Pauli(orbTempZ[i], orbTempZ[i], 1, -1.0, false); //runs
            MSG_INFO("spinorbit before pauli");
            mrcpp::apply_Pauli(orbTempX2, orbTempX[i], 1, -1.0, false);//test debug test
            mrcpp::apply_Pauli(orbTempY2, orbTempY[i], 2, -1.0, false);//test debug test
            mrcpp::apply_Pauli(orbTempZ2, orbTempZ[i], 3, -1.0, false);//test debug test
            MSG_INFO("spinorbit applied pauli");
            ComplexDouble cmplx_i = {0.0, 1.0};
            termSO[i].add(cmplx_i, orbTempX2);//test debug test
            termSO[i].add(cmplx_i, orbTempY2);//test debug test
            termSO[i].add(cmplx_i, orbTempZ2);//test debug test
            // Orbital termSOtemp;                                                                                                                                                                                                        
            // mrcpp::add(termSOtemp, cmplx_i, orbTempX2, cmplx_i, orbTempY2, -1.0, false);                                                                                                                                               
            // mrcpp::add(termSO[i], 1.0, termSOtemp, cmplx_i, orbTempZ2, -1.0, false);  

            //Multiplying by the prefactors. The factor 1/2 comes from the definition of chi, whereas 1/2c^2 comes from the elimination of the small component
            termSO[i].rescale(1 / (2*two_cc)); //could be replaced by simply changing the factor in the addition to the argument

            MSG_INFO("spinorbit orb="<<i << " " << orbTempX2.getSquareNorm()<< " " << orbTempY2.getSquareNorm()<< " " << orbTempZ2.getSquareNorm());
            termSO[i].calcSquareNorm();
            MSG_INFO("SO strength="<< termSO[i].getSquareNorm() << " expectation value for orbital "<< i << ": "<< dot(Phi[i], termSO[i]));
        }
        operSOX.clear();
        operSOY.clear();
        operSOZ.clear();
    }
    // MSG_INFO("Spinorb done");

    //What is this useful for? We don't use them at all.
    auto normsOne = orbital::get_norms(termOne);
    auto normsTwo = orbital::get_norms(termTwo);
    auto normsThree = orbital::get_norms(termThree);
    // auto normsSO = (Phi[0].Ncomp() > 1 ? orbital::get_norms(termSO) : DoubleVector(termSO.size(), 0.0)); //Does not compile. Need to check the constructor for the DoubleVector
    auto normsPsi = orbital::get_norms(Psi);

    // Add up all the terms
    Timer t_add;
    OrbitalVector arg = orbital::deep_copy(termOne);
    for (int i = 0; i < arg.size(); i++) {
        if (not mrcpp::mpi::my_func(arg[i])) continue;
        arg[i].add(1.0, termTwo[i]);
        arg[i].add(1.0, termThree[i]);
        MSG_INFO("pre spinorb: arg["<< i << "] is complex="<<arg[i].iscomplex() << " is real="<<arg[i].isreal());
        if ((Phi[0].Ncomp() > 1) and isZora()) arg[i].add(1.0, termSO[i]); //spin-orbit coupling. Is zero for scalar functions.
        MSG_INFO("post spinorb: arg["<< i << "] is complex="<<arg[i].iscomplex() << " is real="<<arg[i].isreal());
        arg[i].add(1.0, Psi[i]);
    }
    mrcpp::print::time(2, "Adding contributions", t_add);

    operThree.clear();
    operOne.clear();

    Timer t_kappa;
    mrchem::OrbitalVector out = chi_m1(arg);
    for (int i = 0; i < arg.size(); i++) {
        if (not mrcpp::mpi::my_func(out[i])) continue;
        out[i].add(1.0, arg[i]);
    }
    mrcpp::print::time(2, "Applying kappa inverse", t_kappa);
    return out;
}

// Non-relativistic Helmholtz argument
OrbitalVector FockBuilder::buildHelmholtzArgumentNREL(OrbitalVector &Phi, OrbitalVector &Psi) {
    // Get necessary operators
    RankZeroOperator &V = this->potential();

    // Compute OrbitalVectors
    Timer t_pot;
    OrbitalVector termOne = V(Phi);

    mrcpp::print::time(2, "Computing potential term", t_pot);

    // Add up all the terms
    Timer t_add;
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mrcpp::mpi::my_func(out[i])) continue;
        out[i].add(1.0, Psi[i]);
    };
    mrcpp::print::time(2, "Adding contributions", t_add);
    return out;
}

void FockBuilder::setZoraType(bool has_nuc, bool has_coul, bool has_xc, bool is_azora) {
    this->zora_has_nuc = has_nuc;
    this->zora_has_coul = has_coul;
    this->zora_has_xc = has_xc;
    this->zora_is_azora = is_azora;
}

std::shared_ptr<QMPotential> FockBuilder::collectZoraBasePotential() {
    Timer timer;
    auto vz = std::make_shared<QMPotential>(1, false); // normal way
    // auto vz = std::make_shared<QMPotential>(1, false, 2); //TODO: Inclure la quantité de composantes des orbitales dans ce constructeur //Probablement un problème dans le constructeur ici dans le cas où on a 2 composantes? 
    vz->alloc(1, true);
    if (zora_has_nuc) {
        // MSG_INFO("Collecting nuclear potential for ZORA base potential");
        if (getNuclearOperator() != nullptr) {
            // MSG_INFO("Collecting nuclear potential for ZORA base potential tut ");
            auto &vnuc = static_cast<QMPotential &>(getNuclearOperator()->getRaw(0, 0));
            if (not vnuc.hasReal()) MSG_ERROR("ZORA: Adding empty nuclear potential");
            // MSG_INFO("Adding nuclear potential to ZORA base potential");
            vz->add(1.0, vnuc);
        } else {
            MSG_ERROR("ZORA: Nuclear requested but not available");
        }
    }
    if (zora_has_coul) {
        // MSG_INFO("Collecting Coulomb potential for ZORA base potential");
        if (getCoulombOperator() != nullptr) {
            auto &coul = static_cast<QMPotential &>(getCoulombOperator()->getRaw(0, 0));
            // if (not coul.hasReal()) MSG_INFO("ZORA: Adding empty Coulomb potential");
            vz->add(1.0, coul);
        } else {
            MSG_ERROR("ZORA: Coulomb requested but not available");
        }
    }
    if (zora_has_xc) {
        // MSG_INFO("Collecting XC potential for ZORA base potential");
        if (getXCOperator() != nullptr) {
            getXCOperator()->setSpin(SPIN::Paired);
            auto &xc = static_cast<QMPotential &>(getXCOperator()->getRaw(0, 0));
            if (not xc.hasReal()) MSG_ERROR("ZORA: Adding empty XC potential");
            vz->add(1.0, xc);
            getXCOperator()->clearSpin();
        } else {
            MSG_ERROR("ZORA: XC requested but not available");
        }
    }
    print_utils::qmfunction(2, "ZORA operator (base)", *vz, timer);
    return vz;
}

} // namespace mrchem
