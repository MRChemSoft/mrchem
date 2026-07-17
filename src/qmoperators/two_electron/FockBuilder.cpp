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

#include "pseudopotential/projectorOperator.h"

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
    if (this->nuc != nullptr) this->V += (*this->nuc);
    if (this->coul != nullptr) this->V += (*this->coul); 
    if (this->ex != nullptr) this->V -= this->exact_exchange * (*this->ex);
    if (this->xc != nullptr) this->V += (*this->xc);
    if (this->ext != nullptr) this->V += (*this->ext);
    if (this->Ro != nullptr) this->V -= (*this->Ro);
    if (this->pp_projector != nullptr) this->V += (*this->pp_projector);
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



    auto plevel = Printer::getPrintLevel();
    if (plevel >= 2) {
        mrcpp::print::header(2, "Building Fock operator");
        mrcpp::print::value(2, "Precision", prec, "(rel)", 5);
        mrcpp::print::separator(2, '-');
    }

    this->prec = prec;
    if (this->mom != nullptr) this->momentum().setup(prec);
    this->potential().setup(prec);
    this->perturbation().setup(prec); //TODO: uncomment when response is implemented for multiple components


    if (isZora()) {
        Timer t_zora; 
        double c = getLightSpeed();
        mrcpp::print::header(3, "Building ZORA operators");
        mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
        mrcpp::print::value(3, "Light speed", c, "(au)", 5);
        mrcpp::print::separator(3, '-');
        auto vz = collectZoraBasePotential();
        // chi = kappa - 1. See ZoraOperator.h for more information.
        this->chi = std::make_shared<ZoraOperator>(*vz, c, prec, false);
        this->chi_inv = std::make_shared<ZoraOperator>(*vz, c, prec, true);
        this->zora_base = RankZeroOperator(vz);
        this->chi->setup(prec);
        this->chi_inv->setup(prec);
        this->zora_base.setup(prec);
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
    double E_nl = 0.0;   // Non-local pseudopotential energy

    // Nuclear part
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
    if (isZora() || isAZora()) {
        // we assume that at least one orbital is owned by this MPI (TODO: allreduce)
        int Ncomponents = 1;
        for (int i = 0; i < Phi.size(); i++) {
            if (!mrcpp::mpi::my_func(i)) continue;
            Ncomponents = std::max(Ncomponents, Phi[i].Ncomp()); //assumes all owned orbitals have the same number of components
        }
        Ncomponents = mrcpp::mpi::allreduce_max(Ncomponents, mrcpp::mpi::comm_wrk);
        bool spinorial = ( Ncomponents > 1);
        //second term doesn't inclue Pauli matrices (i.e. spinorial is false) because (σ·p)(σ·p) = p^2
        E_kin = qmoperator::calc_kinetic_trace(momentum(), *this->chi, Phi, spinorial).real() + qmoperator::calc_kinetic_trace(momentum(), Phi);
    } else {
        E_kin = qmoperator::calc_kinetic_trace(momentum(), Phi);
    }

    // Electronic part
    if (this->nuc != nullptr) { E_en = this->nuc->trace(Phi).real(); }

    if (this->coul != nullptr) E_ee = 0.5 * this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->ext != nullptr) E_eext = this->ext->trace(Phi).real();
    if (getProjectorOperator() != nullptr) {
        E_nl = this->pp_projector->trace(Phi).real();
    }

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing molecular energy", t_tot);

    return SCFEnergy{E_kin, E_nn, E_en, E_ee, E_x, E_xc, E_next, E_eext, Er_tot, Er_nuc, Er_el, E_nl};
}

/** @brief Compute the Fock matrix F_ij = <bra_i|F|ket_j>
 * 
 */
ComplexMatrix FockBuilder::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix");

    ComplexMatrix T_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    if (isZora() || isAZora()) {
        //If we have spinors, the kinetic operator is of the form (σ·p)V(σ·p), with σ being a Pauli matrix.
        //What this boolean does is enabling the application of the Pauli matrices along the x,y,z momentum operators.
        //NOTE! The second term does not change from being spinorial; (σ·p)(σ·p) = p^2 using the Dirac identity.
        int Ncomponents = 0;
        for (int i = 0; i < bra.size(); i++) {
            if (!mrcpp::mpi::my_func(i)) continue;
            Ncomponents = std::max(Ncomponents, bra[i].Ncomp());
        }
        for (int i = 0; i < ket.size(); i++) {
            if (!mrcpp::mpi::my_func(i)) continue;
            Ncomponents = std::max(Ncomponents, ket[i].Ncomp());
        }
        Ncomponents = mrcpp::mpi::allreduce_max(Ncomponents, mrcpp::mpi::comm_wrk);
        bool spinorial = (Ncomponents > 1); //assumes all orbitals have the same number of components
        T_mat = qmoperator::calc_kinetic_matrix(momentum(), *this->chi, bra, ket, spinorial) + qmoperator::calc_kinetic_matrix(momentum(), bra, ket);
    } else {
        T_mat = qmoperator::calc_kinetic_matrix(momentum(), bra, ket);
    }

    ComplexMatrix V_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    V_mat += potential()(bra, ket);

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Fock matrix", t_tot);
    return T_mat + V_mat;
}

//Deprecated / dead code, should be removed
ComplexMatrix FockBuilder::kineticMatrix(OrbitalVector &bra, OrbitalVector &ket) {
    if (isZora() || isAZora()) {
        return qmoperator::calc_kinetic_matrix(momentum(), *this->chi, bra, ket) + qmoperator::calc_kinetic_matrix(momentum(), bra, ket);
    } else {
        return qmoperator::calc_kinetic_matrix(momentum(), bra, ket);
    }
}

ComplexMatrix FockBuilder::potentialMatrix(OrbitalVector &bra, OrbitalVector &ket) {
    return potential()(bra, ket);
}

OrbitalVector FockBuilder::buildHelmholtzArgument(double prec, OrbitalVector Phi, ComplexMatrix F_mat, ComplexMatrix L_mat) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Helmholtz argument");

    Timer t_rot;

    OrbitalVector Psi = orbital::rotate(Phi, L_mat - F_mat, prec);
    mrcpp::print::time(2, "Rotating orbitals", t_rot);

    OrbitalVector out;
    if (isZora()) { //currently no support for AZORA in the new argument formulation
        int Ncomponents = 1;
        for (int i = 0; i < Phi.size(); i++) {
            if (!mrcpp::mpi::my_func(i)) continue;
            Ncomponents = std::max(Ncomponents, Phi[i].Ncomp());
        }
        Ncomponents = mrcpp::mpi::allreduce_max(Ncomponents, mrcpp::mpi::comm_wrk);
        //note: the new formulation is practical for either NR or 2C, but 1C requires a different correction than 2C which is not implemented yet
        if (Ncomponents==1) out = buildHelmholtzArgumentZORA(Phi, Psi, F_mat.real().diagonal(), prec); 
        if (Ncomponents>1) out = buildHelmholtzArgumentCompact(Phi, Psi); //test debug alt propagator
    } else if (isAZora()){
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
    MSG_INFO("start c="<< c);

    std::shared_ptr<RankZeroOperator> operThreePtr = nullptr;

    if (isZora()) {
        RankZeroOperator &V_zora = this->zora_base;
        operThreePtr = std::make_shared<RankZeroOperator>(V_zora * chi + V_zora); //original
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

    // Compute OrbitalVectors
    Timer t_1;
    OrbitalVector termOne = operOne(Phi);

    mrcpp::print::time(2, "Computing gradient term", t_1);

    Timer t_2;
    OrbitalVector termTwo = V(Phi);

    mrcpp::print::time(2, "Computing potential term", t_2);

    // Compute transformed orbitals scaled by diagonal Fock elements
    Timer t_3;
    OrbitalVector epsPhi = orbital::deep_copy(Phi);
    for (int i = 0; i < epsPhi.size(); i++) {
        if (not mrcpp::mpi::my_func(epsPhi[i])) continue;
        epsPhi[i].rescale(eps[i] / two_cc);
    }
    OrbitalVector termThree = operThree(epsPhi);
    mrcpp::print::time(2, "Computing rescaled potential term", t_3);

    //spin orbit coupling term, which would be identically 0 for scalar functions.
    OrbitalVector termSO(Phi.size());
    int Ncomponents = 1;
    for (int i = 0; i < Phi.size(); i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        Ncomponents = std::max(Ncomponents, Phi[i].Ncomp());
    }
    if ((Ncomponents == 2) and isZora()) {

        //Manually implementing the cross product appearing in the spin-orbit term
        //NOTE! The multiplication by the Pauli matrices will need to be handled later,
        //      during the application of the cross-product to the orbitals.
        RankZeroOperator operSOX = 0.5 * p(chi)[1]*p[2] - 0.5 * p(chi)[2]*p[1];
        operSOX.setup(prec);
        RankZeroOperator operSOY = 0.5 * p(chi)[2]*p[0] - 0.5 * p(chi)[0]*p[2];
        operSOY.setup(prec);
        RankZeroOperator operSOZ = 0.5 * p(chi)[0]*p[1] - 0.5 * p(chi)[1]*p[0];
        operSOZ.setup(prec);
        //Applying the curls to temporary copies of the orbitals 
        //NOTE! Extremely inefficient!
        OrbitalVector orbTempX = operSOX(Phi);
        OrbitalVector orbTempY = operSOY(Phi);
        OrbitalVector orbTempZ = operSOZ(Phi);
        //adding the contributions together. Note that the coefficient is purely imaginary.
        for (int i = 0; i < Phi.size(); i++) {
            mrcpp::apply_gamma(orbTempX[i], 1);
            mrcpp::apply_gamma(orbTempY[i], 2);
            mrcpp::apply_gamma(orbTempZ[i], 3);

            ComplexDouble cmplx_i = {0.0, 1.0}; 
            termSO[i].add(cmplx_i, orbTempX[i]);
            termSO[i].add(cmplx_i, orbTempY[i]);
            termSO[i].add(cmplx_i, orbTempZ[i]);

            termSO[i].crop(prec);
        }
        operSOX.clear();
        operSOY.clear();
        operSOZ.clear();
    }

    //What is this useful for? We don't use them at all.
    auto normsOne = orbital::get_norms(termOne);
    auto normsTwo = orbital::get_norms(termTwo);
    auto normsThree = orbital::get_norms(termThree);
    // auto normsSO = (Phi[0].Ncomp() > 1 ? orbital::get_norms(termSO) : DoubleVector(termSO.size(), 0.0)); //Don't know why we compute them norms, we don't use them later
    auto normsPsi = orbital::get_norms(Psi);

    // Add up all the terms
    Timer t_add;
    OrbitalVector arg = orbital::deep_copy(termOne);
    for (int i = 0; i < arg.size(); i++) {
        if (not mrcpp::mpi::my_func(arg[i])) continue;
        arg[i].add(1.0, termTwo[i]);
        arg[i].add(1.0, termThree[i]);
        if ((Ncomponents > 1) and isZora()) arg[i].add(1.0, termSO[i]); //spin-orbit coupling. Is zero for scalar functions.
        arg[i].add(1.0, Psi[i]);
    }
    mrcpp::print::time(2, "Adding contributions", t_add);

    operThree.clear();
    operOne.clear();

    Timer t_kappa;
    mrchem::OrbitalVector out = chi_m1(arg);
    //chi^{-1} is defined as (kappa^{-1})-1 so we need to add the argument back
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

// Propagator integrand construction (inhomogeneous part of the Helmholtz equation)
OrbitalVector FockBuilder::buildHelmholtzArgumentCompact(OrbitalVector &Phi, OrbitalVector &Psi) {
    // Get necessary operators
    double c = getLightSpeed();
    double two_cc = 2.0 * c * c;
    MomentumOperator &p = momentum();
    RankZeroOperator &V = potential();
    RankZeroOperator &chi = *this->chi;
    // RankZeroOperator &chi_m1 = *this->chi_inv; //not needed in this formulation

    // Compute OrbitalVectorso
    Timer t_pot;
    OrbitalVector termOne = V(Phi);

    OrbitalVector termTwo(Phi.size());//compute ZORA correction term (σ·p)χ(σ·p)|ψ>
    //Need to know if we have to apply the σ matrices in a MPI-safe way
    int Ncomponents = 1;
    for (int i = 0; i < Phi.size(); i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        Ncomponents = std::max(Ncomponents, Phi[i].Ncomp());
    }
    Ncomponents = mrcpp::mpi::allreduce_max(Ncomponents, mrcpp::mpi::comm_wrk);
    // boolean to toggle the application of the Pauli/gamma matrices to the spinors
    bool rotate_spin = (Ncomponents > 1); 
    for (int i = 0; i < Phi.size(); i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        // apply (σ_d p_d) to |ψ> (no summation yet)
        std::vector<Orbital> nabla_Phi = p(Phi[i], rotate_spin);
        //sum it up in one temporary term
        Orbital termTwo_tmp; //holds (σ·p)|ψ>
        termTwo_tmp.add({1.0, 0.0}, nabla_Phi[0]);
        termTwo_tmp.add({1.0, 0.0}, nabla_Phi[1]);
        termTwo_tmp.add({1.0, 0.0}, nabla_Phi[2]);
        
        // Free memory space by discarding no longer relevant trees. Should help mitigate the memory usage spike from this function
        for (int dim=0; dim<3; dim++) nabla_Phi[dim].free();

        // apply χ to (σ·p)|ψ>
        termTwo_tmp = chi(termTwo_tmp); 
        
        // apply (σ_d p_d) to χ(σ·p)|ψ> (no summation yet)
        std::vector<Orbital> nabla_chi_nabla_Phi = p(termTwo_tmp, rotate_spin);
        // sum it up in the final term
        termTwo[i].add({1.0, 0.0}, nabla_chi_nabla_Phi[0]);
        termTwo[i].add({1.0, 0.0}, nabla_chi_nabla_Phi[1]);
        termTwo[i].add({1.0, 0.0}, nabla_chi_nabla_Phi[2]);
    }

    mrcpp::print::time(2, "Computing potential term", t_pot);

    // Add up all the terms to form the inhomogeneous part of the Helmholtz equation
    Timer t_add;
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mrcpp::mpi::my_func(out[i])) continue;
        out[i].add(0.5, termTwo[i]); 
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
    vz->alloc(1, true);
    if (zora_has_nuc) {
        if (getNuclearOperator() != nullptr) {
            auto &vnuc = static_cast<QMPotential &>(getNuclearOperator()->getRaw(0, 0));
            if (not vnuc.hasReal()) MSG_ERROR("ZORA: Adding empty nuclear potential");
            vz->add(1.0, vnuc);
        } else {
            MSG_ERROR("ZORA: Nuclear requested but not available");
        }
    }
    if (zora_has_coul) {
        if (getCoulombOperator() != nullptr) {
            auto &coul = static_cast<QMPotential &>(getCoulombOperator()->getRaw(0, 0));
            if (not coul.hasReal()) MSG_INFO("ZORA: Adding empty Coulomb potential");
            vz->add(1.0, coul);
        } else {
            MSG_ERROR("ZORA: Coulomb requested but not available");
        }
    }
    if (zora_has_xc) {
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
