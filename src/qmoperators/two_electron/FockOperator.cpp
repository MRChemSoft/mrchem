#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "FockOperator.h"
#include "KineticOperator.h"
#include "NuclearOperator.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "ElectricFieldOperator.h"
#include "XCFunctional.h"
#include "XCOperator.h"
#include "SCFEnergy.h"
#include "utils/math_utils.h"
#include "chemistry/chem_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

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
FockOperator::FockOperator(KineticOperator  *t,
                           NuclearOperator  *v,
                           CoulombOperator  *j,
                           ExchangeOperator *k,
                           XCOperator       *xc,
                           ElectricFieldOperator *ext)
        : kin(t),
          nuc(v),
          coul(j),
          ex(k),
          xc(xc),
          ext(ext){
}

/** @brief build the Fock operator once all contributions are in place
 * 
 */
void FockOperator::build() {
    if (this->kin  != 0) this->T  = *this->kin;
    if (this->nuc  != 0) this->V  = *this->nuc;
    if (this->coul != 0) this->V += *this->coul;
    if (this->ex   != 0) this->V -= *this->ex;
    if (this->xc   != 0) this->V += *this->xc;
    if (this->ext  != 0) this->V += *this->ext;
    RankZeroTensorOperator &F = (*this);
    F = this->T + this->V;
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will call the setup function of all underlying operators, and in particular
 * it will compute the internal exchange if there is an ExchangeOperator.
 */
void FockOperator::setup(double prec) {
    Timer timer;
    Printer::printHeader(0, "Setting up Fock operator");
    Printer::printDouble(0, "Precision", prec, 5);
    Printer::printSeparator(0, '-');
    this->T.setup(prec);
    this->V.setup(prec);
    if (this->ex != 0) this->ex->setupInternal(prec);
    timer.stop();
    Printer::printFooter(0, timer, 2);
}

/** @brief clear operator after application
 *
 * This will call the clear function of all underlying operators, and bring them back
 * to the state after construction. The operator can now be reused after another setup.
 */
void FockOperator::clear() {
    this->T.clear();
    this->V.clear();
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
    if (this->ex != 0) this->ex->rotate(U);
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
SCFEnergy FockOperator::trace(OrbitalVector &Phi, const ComplexMatrix &F) {
    double E_nuc = 0.0; // Nuclear repulsion
    double E_el  = 0.0; // Electronic energy
    double E_orb = 0.0; // Orbital energy
    double E_kin = 0.0; // Kinetic energy
    double E_en  = 0.0; // Nuclear-electronic interaction
    double E_ee  = 0.0; // Electronic repulsion
    double E_x   = 0.0; // Exact Exchange
    double E_xc  = 0.0; // Exchange and Correlation
    double E_xc2 = 0.0; // Trace of the XC operator
    double E_ext = 0.0; // External field contribution to the electronic energy
    double E_nex = 0.0; // External field contribution to the nuclear energy

    // Nuclear part
    if (this->nuc != 0) {
        Nuclei &nucs = this->nuc->getNuclei();
        E_nuc = compute_nuclear_repulsion(nucs);
        if (this->ext  != 0) {
            E_nex = this->ext->trace(nucs).real();
            E_nuc += E_nex;
        }
    }

    // Orbital energies
    for (int i = 0; i < Phi.size(); i++) {
        double occ = (double) Phi[i].occ();
        E_orb += occ*F(i,i).real();
    }

    // Electronic part
    if (this->nuc  != 0) E_en  =  this->nuc->trace(Phi).real();
    if (this->coul != 0) E_ee  =  this->coul->trace(Phi).real();
    if (this->ex   != 0) E_x   = -this->ex->trace(Phi).real();
    if (this->xc   != 0) E_xc  =  this->xc->getEnergy();
    if (this->xc   != 0) E_xc2 =  this->xc->trace(Phi).real();
    if (this->ext  != 0) E_ext =  this->ext->trace(Phi).real();

    double E_eex    = E_ee  + E_x;
    double E_orbxc2 = E_orb - E_xc2;
    E_kin = E_orbxc2 - 2.0*E_eex - E_en - E_ext;
    E_el  = E_orbxc2 -     E_eex + E_xc;

    return SCFEnergy(E_nuc, E_el, E_orb, E_kin, E_en, E_ee, E_xc, E_x,
                     E_nex, E_ext);
}

} //namespace mrchem
