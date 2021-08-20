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
#include "FarFieldOperator.h"
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
#include "qmoperators/one_electron/SmearedNuclearOperator.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using KineticOperator_p = std::shared_ptr<mrchem::KineticOperator>;
using NuclearOperator_p = std::shared_ptr<mrchem::NuclearOperator>;
using SmearedNuclearOperator_p = std::shared_ptr<mrchem::SmearedNuclearOperator>;
using CoulombOperator_p = std::shared_ptr<mrchem::CoulombOperator>;
using FarFieldOperator_p = std::shared_ptr<mrchem::FarFieldOperator>;
using ExchangeOperator_p = std::shared_ptr<mrchem::ExchangeOperator>;
using XCOperator_p = std::shared_ptr<mrchem::XCOperator>;
using ElectricFieldOperator_p = std::shared_ptr<mrchem::ElectricFieldOperator>;
using ReactionOperator_p = std::shared_ptr<mrchem::ReactionOperator>;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param t:   Kinetic operator
 * @param v:   Nuclear potential operator
 * @param j:   Coulomb operator
 * @param k:   HF exchange operator
 * @param xc:  Exchange-Correlation operator
 * @param ext: External Field operator
 *
 * Each of the arguments can be NULL, so this operators includes both core Hamiltonian,
 * the Hartree(-Fock) method and (pure/hybrid) Density Functional Theory.
 */
FockOperator::FockOperator(KineticOperator_p t,
                           NuclearOperator_p v,
                           SmearedNuclearOperator_p sv,
                           CoulombOperator_p j,
                           FarFieldOperator_p f,
                           ExchangeOperator_p k,
                           XCOperator_p xc,
                           ElectricFieldOperator_p ext,
                           ReactionOperator_p reo)
        : kin(t)
        , nuc(v)
        , coul(j)
        , ff(f)
        , ex(k)
        , xc(xc)
        , ext(ext)
        , Ro(reo) {}

/** @brief build the Fock operator once all contributions are in place
 *
 */
void FockOperator::build(double exx) {
    this->exact_exchange = exx;
    this->T = RankZeroOperator();
    if (this->kin != nullptr) this->T += (*this->kin);

    this->V = RankZeroOperator();
    if (this->nuc != nullptr) this->V += (*this->nuc);
    if (this->snuc != nullptr) this->V += (*this->snuc);
    if (this->coul != nullptr) this->V += (*this->coul);
    if (this->ff != nullptr) this->V += (*this->ff);
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
    double E_eeff = 0.0; // Far field contribution to the electronic energy
    double E_nff = 0.0;  // Far field contribution to the nuclear energy

    //// PERIODIC ENERGIES ////
    double E_se = 0.0;
    double E_ncor = 0.0;

    // Nuclear part
    if (this->nuc != nullptr) E_nn = chemistry::compute_nuclear_repulsion(nucs);
    if (this->ext != nullptr) E_next = -this->ext->trace(nucs).real();
    if (this->ff != nullptr) E_nff = -0.5 * this->ff->trace(nucs).real();

    // Reaction potential part
    if (this->Ro != nullptr) {
        Er_nuc = 0.5 * this->Ro->getNuclearEnergy();
        Er_tot = 0.5 * this->Ro->getTotalEnergy();
        Er_el = 0.5 * this->Ro->getElectronicEnergy();
    }
    // Electronic part
    if (this->kin != nullptr) E_kin = this->kin->trace(Phi).real();
    if (this->nuc != nullptr) E_en = this->nuc->trace(Phi).real();
    if (this->snuc != nullptr) E_en = this->snuc->trace(Phi).real();
    if (this->coul != nullptr) E_ee = 0.5 * this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->ext != nullptr) E_eext = this->ext->trace(Phi).real();
    if (this->ff != nullptr) E_eeff = this->ff->trace(Phi).real();

    if (this->coul->getPotential()->getRc() > 0.0) {
        E_nn = 0.0;
        E_se = -0.5 * mrcpp::dot<3>(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
        auto rc = this->coul->getPotential()->getRc();
        if (rc <= 0.0) MSG_ABORT("RC has to be positive");
        auto Ig = 10976.0 / (17875.0 * rc);
        for (auto &nuc : nucs) E_ncor += 0.5 * nuc.getCharge() * nuc.getCharge() * (Ig - 12.0 / (5.0 * rc));
    }

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing molecular energy", t_tot);

    return SCFEnergy{E_kin, E_nn, E_en, E_ee, E_x, E_xc, E_next, E_eext, Er_tot, Er_nuc, Er_el, E_nff, E_eeff, E_se, E_ncor};
}

ComplexMatrix FockOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix");

    auto t = this->getKineticOperator();
    auto v = this->potential();

    ComplexMatrix T = ComplexMatrix::Zero(bra.size(), ket.size());
    if (t != nullptr) T += qmoperator::calc_kinetic_matrix(*t, bra, ket);

    ComplexMatrix V = ComplexMatrix::Zero(bra.size(), ket.size());
    if (v.size() > 0) V += v(bra, ket);

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Fock matrix", t_tot);
    return T + V;
}

ComplexMatrix FockOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
