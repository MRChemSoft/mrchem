#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MathUtils.h"

#include "FockOperator.h"
#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "SCFEnergy.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

FockOperator::FockOperator(KineticOperator  *t,
                           NuclearPotential *v,
                           CoulombOperator  *j,
                           ExchangeOperator *k,
                           XCOperator       *xc)
        : kin(t),
          nuc(v),
          coul(j),
          ex(k),
          xc(xc) {
    if (this->kin  != 0) this->T += *this->kin;
    if (this->nuc  != 0) this->V += *this->nuc;
    if (this->coul != 0) this->V += *this->coul;
    if (this->ex   != 0) this->V -= *this->ex;
//    if (this->xc   != 0) this->V += *this->xc;

    RankZeroTensorOperator &F = (*this);
    F = this->T + this->V;
}

void FockOperator::setup(double prec) {
    Timer timer;
    Printer::printHeader(0, "Setting up Fock operator");
    Printer::printDouble(0, "Precision", prec);
    Printer::printSeparator(0, '-');
    this->T.setup(prec);
    this->V.setup(prec);
    if (this->ex != 0) this->ex->setupInternal(prec);
    timer.stop();
    Printer::printFooter(0, timer, 2);
}

void FockOperator::clear() {
    this->T.clear();
    this->V.clear();
}

void FockOperator::rotate(const ComplexMatrix &U) {
    if (this->ex != 0) this->ex->rotate(U);
}

SCFEnergy FockOperator::trace(OrbitalVector &Phi, const ComplexMatrix &F) {
    double E_nuc = 0.0;
    double E_el  = 0.0;
    double E_orb = 0.0;
    double E_kin = 0.0;
    double E_en  = 0.0;
    double E_ee  = 0.0;
    double E_x   = 0.0;
    double E_xc  = 0.0;
    double E_xc2 = 0.0;

    // Nuclear part
    if (this->nuc != 0) {
        Nuclei &nucs = this->nuc->getNuclei();
        int nNucs = nucs.size();
        for (int i = 0; i < nNucs; i++) {
            for (int j = i+1; j < nNucs; j++) {
                const Nucleus &nuc_i = nucs[i];
                const Nucleus &nuc_j = nucs[j];
                double Z_i = nuc_i.getCharge();
                double Z_j = nuc_j.getCharge();
                const double *R_i = nuc_i.getCoord();
                const double *R_j = nuc_j.getCoord();
                double R_ij = mrcpp::MathUtils::calcDistance(3, R_i, R_j);
                E_nuc += (Z_i*Z_j)/R_ij;
            }
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
    //if (this->xc   != 0) E_xc  =  this->xc->getEnergy();
    //if (this->xc   != 0) E_xc2 =  this->xc->trace(Phi).real();

    double E_eex    = E_ee  + E_x;
    double E_orbxc2 = E_orb - E_xc2;
    E_kin = E_orbxc2 - 2.0*E_eex - E_en;
    E_el  = E_orbxc2 - E_eex     + E_xc;

    return SCFEnergy(E_nuc, E_el, E_orb, E_kin, E_en, E_ee, E_xc, E_x);
}

} //namespace mrchem
