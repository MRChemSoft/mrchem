#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "utils/math_utils.h"

#include "NuclearOperator.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

NuclearPotential::NuclearPotential(const Nuclei &nucs, double prec)
        : QMPotential(1, mpi::share_nuc_pot) {
    int oldprec = Printer::setPrecision(5);
    Printer::printHeader(0, "Setting up nuclear potential");
    println(0, " Nr  Element         Charge        Precision     Smoothing ");
    Printer::printSeparator(0, '-');

    double c = 0.00435 * prec;
    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        double Z = nuc.getCharge();
        double Z_5 = pow(Z, 5.0);
        double smooth = pow(c / Z_5, 1.0 / 3.0);
        this->func.push_back(nuc, smooth);

        std::stringstream symbol;
        symbol << nuc.getElement().getSymbol();
        symbol << "  ";
        printout(0, std::setw(3) << i + 1 << "     ");
        printout(0, symbol.str()[0] << symbol.str()[1]);
        printout(0, std::setw(22) << Z);
        printout(0, std::setw(14) << prec);
        printout(0, std::setw(14) << smooth << std::endl);
    }
    Printer::printSeparator(0, '=', 2);
    Printer::setPrecision(oldprec);
}

void NuclearPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    QMPotential &V = *this;

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    qmfunction::project(V, this->func, NUMBER::Real, this->apply_prec);
    timer.stop();
    int n = V.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(0, "Nuclear potential", n, t);
}

void NuclearPotential::clear() {
    free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();    // apply_prec = -1
}

/** @brief computes the interaction energy of the nuclear potential with a second set of nuclei.
 * 
 * @param[in] nucs the set of nuclei
 *
 * Note: this function is not suited to compute the nuclear self-energy
 *
 */
double NuclearOperator::trace(const Nuclei &nucs) {
    MSG_WARN("This routine has never been tested!");
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int i = 0; i < nNucs; i++) {
        const Nucleus &nuc_i = nucs[i];
        double Z_i = nuc_i.getCharge();
        const mrcpp::Coord<3> &R_i = nuc_i.getCoord();
        E_nuc += Z_i * this->r_m1.evalf(R_i);
    }
    return E_nuc;
}

} //namespace mrchem
