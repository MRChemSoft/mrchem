#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NuclearOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

NuclearPotential::NuclearPotential(const Nuclei &nucs, double prec)
        : QMPotential(1) {
    int oldprec = Printer::setPrecision(5);
    Printer::printHeader(0, "Setting up nuclear potential");
    println(0, " Nr  Element         Charge        Precision     Smoothing ");
    Printer::printSeparator(0, '-');

    double c = 0.00435*prec;
    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        double Z = nuc.getCharge();
        double Z_5 = pow(Z, 5.0);
        double smooth = pow(c/Z_5, 1.0/3.0);
        this->func.push_back(nuc, smooth);

        std::stringstream symbol;
        symbol << nuc.getElement().getSymbol();
        symbol << "  ";
        printout(0, std::setw(3) << i+1 << "     ");
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

    if (hasReal()) MSG_ERROR("Potential not properly cleared");
    if (hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    alloc(NUMBER::Real);
    mrcpp::build_grid(this->real(), this->func);
    mrcpp::project(this->apply_prec, this->real(), this->func);
    timer.stop();

    int n = getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Nuclear potential", n, t);
}

void NuclearPotential::clear() {
    free();           // delete FunctionTree pointers
    clearApplyPrec(); // apply_prec = -1
}

    /** @brief computes the repulsion self energy of a set of nuclei
     * 
     * @param[in] nucs the set of nuclei
     *
     */
    
double NuclearOperator::trace(const Nuclei &nucs) {
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int i = 0; i < nNucs; i++) {
        for (int j = i+1; j < nNucs; j++) {
            const Nucleus &nuc_i = nucs[i];
            const Nucleus &nuc_j = nucs[j];
            double Z_i = nuc_i.getCharge();
            double Z_j = nuc_j.getCharge();
            const double *R_i = nuc_i.getCoord();
                const double *R_j = nuc_j.getCoord();
                double R_ij = math_utils::calc_distance(R_i, R_j);
                E_nuc += (Z_i*Z_j)/R_ij;
        }
    }
    return E_nuc;
}

} //namespace mrchem
