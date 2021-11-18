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

    if (this->sqrt_vz != nullptr) this->sqrt_zora += (*this->sqrt_vz);
    if (this->mod_vz != nullptr) this->mod_zora += (*this->mod_vz);
    if (this->scaled_vz != nullptr) this->scaled_zora += (*this->scaled_vz);
    if (this->inv_vz != nullptr) this->inv_zora += (*this->inv_vz);
    if (this->one_minus_vz != nullptr) this->one_minus_zora += *(this->one_minus_vz);

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
    this->sqrt_zora_pot().setup(prec);
    this->mod_zora_pot().setup(prec);
    this->scaled_zora_pot().setup(prec);
    this->inv_zora_pot().setup(prec);
    this->one_minus_zora_pot().setup(prec);
    this->zora().setup(prec);
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
        this->sqrt_zora_pot().clear();
        this->mod_zora_pot().clear();
        this->scaled_zora_pot().clear();
        this->inv_zora_pot().clear();
        this->one_minus_zora_pot().clear();
        this->zora().clear();
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
        switch (this->zoraKineticAlgorithm) {
            case 0:  // The "straightforward" way
                {
                T_mat = qmoperator::calc_kinetic_matrix(momentum(), zora(), bra, ket);
                break;
                }
            case 1:  // Trying to exploit symmetry in laplacian operator
                {
                RankZeroOperator sqrtK(this->sqrt_zora_pot());
                T_mat = qmoperator::calc_kinetic_matrix_symmetrized(momentum(), sqrtK, bra, ket);
                break;
                }
            case 2:  // Using transformed orbitals, Florian Bischoff
                {
                RankZeroOperator sqrtK = this->sqrt_zora_pot();
                RankZeroOperator modK = this->mod_zora_pot();

                OrbitalVector kKet = sqrtK(ket);
                OrbitalVector kBra = sqrtK(bra);

                ComplexMatrix T_1 = ComplexMatrix::Zero(kBra.size(), kKet.size());
                ComplexMatrix T_2 = ComplexMatrix::Zero(kBra.size(), kKet.size());
                
                T_1 = qmoperator::calc_kinetic_matrix(momentum(), kBra, kKet);
                T_2 = modK(kBra, kKet);

                T_mat = T_1 + 0.5 * T_2;
                break;
                }
            case 3:  // Using 1-kappa (to be used with Take 4)
                {
                T_mat = qmoperator::calc_kinetic_matrix(momentum(), bra, ket) - qmoperator::calc_kinetic_matrix(momentum(), zora(), bra, ket);
                break;
                }
        }
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

// Take 2 in notes on Overleaf
OrbitalVector FockOperator::buildHelmholtzArgumentTake2(OrbitalVector &Phi, OrbitalVector &Psi, double prec) {
    // Get necessary operators
    RankZeroOperator O;
    RankZeroOperator &V = this->potential();
    ZoraOperator &kappa = this->zora();
    auto &diff = kappa.derivative;
    NablaOperator nabla(diff);
    KineticOperator T(diff);
    IdentityOperator I;
    RankOneOperator<3> dKappa = nabla(kappa);

    // Construct and set up
    O += (kappa - I) * T;
    O += -0.5 * tensor::dot(dKappa, nabla);
    O += V;
    O.setup(prec);

    // Apply on orbitals
    OrbitalVector termOne = O(Phi);

    // Add Psi
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;
        out[i].add(1.0, Psi[i]);
    };

    // Clear operator
    O.clear();
    return out;
}

// Take 3 in notes on Overleaf
OrbitalVector FockOperator::buildHelmholtzArgumentTake3(OrbitalVector &Phi, OrbitalVector &Psi, DoubleVector eps) {
    // Get necessary operators
    RankZeroOperator &V = this->potential();                 // V
    ZoraOperator &kappa = this->zora();                      // kappa
    RankZeroOperator &zfacKappa = this->scaled_zora_pot();   // Vz / 2c^2
    RankZeroOperator &divKappa = this->mod_zora_pot();       // div(sqrt(kappa)) / sqrt(kappa)
    RankZeroOperator &sqKappa = this->sqrt_zora_pot();       // sqrt(kappa)
    RankZeroOperator &invKappa = this->inv_zora_pot();       // 1 / kappa

    // Compute transformed orbitals scaled by diagonal Fock elements
    OrbitalVector epsPhi = orbital::deep_copy(Phi);
    for (int i = 0; i < epsPhi.size(); i++) {
        if (not mpi::my_orb(epsPhi[i])) continue;
        epsPhi[i].rescale(eps[i]);
    };

    // Compute OrbitalVectors
    OrbitalVector termOne = (0.5 * divKappa)(Phi);
    OrbitalVector termTwo = (V * invKappa)(Phi);
    OrbitalVector termThree = zfacKappa(epsPhi);
    OrbitalVector termFour = invKappa(Psi);

    // Add up all the terms
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;
        out[i].add(1.0, termTwo[i]);
        out[i].add(1.0, termThree[i]);
        out[i].add(1.0, termFour[i]);
    };
    return out;
}

// Take 4 in notes on Overleaf
OrbitalVector FockOperator::buildHelmholtzArgumentTake4(OrbitalVector &Phi, OrbitalVector &Psi, double prec) {
    // Get necessary operators
    RankZeroOperator O;
    RankZeroOperator &V = this->potential();
    RankZeroOperator &omk = this->one_minus_zora_pot();
    auto &diff = this->zora().derivative;
    NablaOperator nabla(diff);
    KineticOperator T(diff);
    RankOneOperator<3> dOmk = nabla(omk);

    // Construct and set up
    O += -1.0 * omk * T;
    O += 0.5 * tensor::dot(dOmk, nabla);
    O += V;
    O.setup(prec);

    // Apply on orbitals
    OrbitalVector termOne = O(Phi);

    // Add Psi
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;
        out[i].add(1.0, Psi[i]);
    };

    // Clear operator
    O.clear();
    return out;
}

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
