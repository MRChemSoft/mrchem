/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "FockOperator.h"
#include "ReactionOperator.h"
#include "XCOperator.h"
#include "chemistry/chemistry_utils.h"
#include "properties/SCFEnergy.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include "qmoperators/one_electron/IdentityOperator.h"
#include "qmoperators/one_electron/ZoraKineticOperator.h"
#include "qmoperators/qmoperator_utils.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief build the Fock operator once all contributions are in place
 *
 */
void FockOperator::build(double exx) {
    this->exact_exchange = exx;

    this->T = RankZeroOperator();
    if (this->mom != nullptr) {
        if (isZora()) {
            this->T += ZoraKineticOperator(momentum(), zora());
        } else {
            this->T += KineticOperator(momentum());
        }
    }

    if (this->vz_over_2cc != nullptr) this->zora_over_2cc = (*this->vz_over_2cc);
    if (this->vz_inv != nullptr) this->zora_inv = (*this->vz_inv);
    if (this->vz_grad != nullptr) this->zora_grad = *(this->vz_grad);

    this->V = RankZeroOperator();
    if (this->nuc != nullptr) this->V += (*this->nuc);
    if (this->coul != nullptr) this->V += (*this->coul);
    if (this->ex != nullptr) this->V -= this->exact_exchange * (*this->ex);
    if (this->xc != nullptr) this->V += (*this->xc);
    if (this->ext != nullptr) this->V += (*this->ext);
    if (this->Ro != nullptr) this->V -= (*this->Ro);

    RankZeroOperator &F = (*this);
    F = this->kinetic() + this->potential();
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will call the setup function of all underlying operators, and in particular
 * it will compute the internal exchange if there is an ExchangeOperator.
 */
void FockOperator::setup(double prec) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Building Fock operator");
    mrcpp::print::value(2, "Precision", prec, "(rel)", 5);
    mrcpp::print::separator(2, '-');
    this->kinetic().setup(prec);
    this->potential().setup(prec);
    this->perturbation().setup(prec);
    
    if (isZora()) {
        this->zora().setup(prec);
        this->zora_pot_over_2cc().setup(prec);
        this->zora_pot_inv().setup(prec);
        this->zora_pot_grad().setup(prec);
    };

    t_tot.stop();
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_tot);
}

/** @brief clear operator after application
 *
 * This will call the clear function of all underlying operators, and bring them back
 * to the state after construction. The operator can now be reused after another setup.
 */
void FockOperator::clear() {
    this->kinetic().clear();
    this->potential().clear();
    this->perturbation().clear();
    
    if (isZora()) {
        this->zora().clear();
        this->zora_pot_over_2cc().clear();
        this->zora_pot_inv().clear();
        this->zora_pot_grad().clear();
    };
}

/** @brief rotate orbitals of two-electron operators
 *
 * @param U: unitary transformation matrix
 *
 * This function should be used in case the orbitals are rotated *after* the FockOperator
 * has been setup. In particular the ExchangeOperator needs to rotate the precomputed
 * internal exchange potentials.
 */
void FockOperator::rotate(const ComplexMatrix &U) {
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
SCFEnergy FockOperator::trace(OrbitalVector &Phi, const Nuclei &nucs) {
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
    if (this->nuc != nullptr) E_nn = chemistry::compute_nuclear_repulsion(nucs);
    if (this->ext != nullptr) E_next = -this->ext->trace(nucs).real();

    // Reaction potential part
    if (this->Ro != nullptr) {
        Er_nuc = 0.5 * this->Ro->getNuclearEnergy();
        Er_tot = 0.5 * this->Ro->getTotalEnergy();
        Er_el = 0.5 * this->Ro->getElectronicEnergy();
    }

    // Kinetic part
    if (isZora()) {
        E_kin = qmoperator::calc_kinetic_trace(momentum(), zora(), Phi).real();
    } else {
        E_kin = qmoperator::calc_kinetic_trace(momentum(), Phi);
    }

    // Electronic part
    if (this->nuc != nullptr) E_en = this->nuc->trace(Phi).real();
    if (this->coul != nullptr) E_ee = 0.5 * this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->ext != nullptr) E_eext = this->ext->trace(Phi).real();
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing molecular energy", t_tot);

    return SCFEnergy{E_kin, E_nn, E_en, E_ee, E_x, E_xc, E_next, E_eext, Er_tot, Er_nuc, Er_el};
}

ComplexMatrix FockOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix");

    ComplexMatrix T_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    if (isZora()) {
        T_mat = qmoperator::calc_kinetic_matrix(momentum(), zora(), bra, ket);
    } else {
        T_mat = qmoperator::calc_kinetic_matrix(momentum(), bra, ket);
    }

    ComplexMatrix V_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    V_mat += potential()(bra, ket);

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Fock matrix", t_tot);
    return T_mat + V_mat;
}

ComplexMatrix FockOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

// Take 1 in notes on Overleaf
OrbitalVector FockOperator::buildHelmholtzArgumentTake1(OrbitalVector &Phi, OrbitalVector &Psi, DoubleVector eps, double prec) {
    // Get necessary operators
    RankZeroOperator &V = this->potential();                   // V
    ZoraOperator &kappa = this->zora();                        // kappa
    RankZeroOperator &divby2cc = this->zora_pot_over_2cc();    // Vz / 2c^2
    RankZeroOperator &invKappa = this->zora_pot_inv();         // 1 / kappa
    RankOneOperator<3> &gradKappa = this->zora_pot_grad();     // gradient of kappa
    NablaOperator nabla(kappa.derivative);
    RankZeroOperator gradKappaGrad = tensor::dot(gradKappa, nabla);
    gradKappaGrad.setup(prec);

    // Compute transformed orbitals scaled by diagonal Fock elements
    OrbitalVector epsPhi = orbital::deep_copy(Phi);
    for (int i = 0; i < epsPhi.size(); i++) {
        if (not mpi::my_orb(epsPhi[i])) continue;
        epsPhi[i].rescale(eps[i]);
    };

    // Compute OrbitalVectors
    OrbitalVector termOne = (-0.5 * gradKappaGrad)(Phi);
    OrbitalVector termTwo = V(Phi);
    OrbitalVector termThree = (kappa * divby2cc)(epsPhi);

    // Add up all the terms
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;
        out[i].add(1.0, termTwo[i]);
        out[i].add(1.0, termThree[i]);
        out[i].add(1.0, Psi[i]);
    };
    
    gradKappaGrad.clear();
    return invKappa(out);
}

// Non-relativistic Helmholtz argument
OrbitalVector FockOperator::buildHelmholtzArgument(OrbitalVector &Phi, OrbitalVector &Psi) {
    // Get necessary operators
    RankZeroOperator &V = this->potential();

    // Compute OrbitalVectors
    OrbitalVector termOne = V(Phi);

    // Add up all the terms
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;
        out[i].add(1.0, Psi[i]);
    };
    return out;
}

} // namespace mrchem
