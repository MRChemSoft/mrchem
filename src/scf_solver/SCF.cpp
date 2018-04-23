#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "SCF.h"
#include "Orbital.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

SCF::SCF(HelmholtzVector &h)
        : maxIter(-1),
          rotation(0),
          canonical(true),
          orbThrs(-1.0),
          propThrs(-1.0),
          orbPrec{-1.0, -1.0, -1.0},
          helmholtz(&h) {
}

void SCF::setThreshold(double orb_thrs, double prop_thrs) {
    this->orbThrs = orb_thrs;
    this->propThrs = prop_thrs;
}

void SCF::setOrbitalPrec(double init, double final) {
    this->orbPrec[0] = init;
    this->orbPrec[1] = init;
    this->orbPrec[2] = final;
}

double SCF::adjustPrecision(double error) {
    if (this->orbPrec[0] > 0.0 ) this->orbPrec[0] *= 0.75;
    this->orbPrec[0] = std::min(10.0*error*error, this->orbPrec[0]);
    this->orbPrec[0] = std::max(this->orbPrec[0], this->orbPrec[2]);

    Printer::printSeparator(0, '=');
    Printer::printDouble(0, "Current precision", this->orbPrec[0], 5);
    Printer::printSeparator(0, '-');
    Printer::printDouble(0, "Orbital threshold", this->orbThrs, 5);
    Printer::printDouble(0, "Property threshold", this->propThrs, 5);
    Printer::printSeparator(0, '=', 2);
    return this->orbPrec[0];
}

void SCF::resetPrecision() {
    this->orbPrec[0] = this->orbPrec[1];
}

bool SCF::checkConvergence(double err_o, double err_p) const {
    bool conv_o = false;
    bool conv_p = false;
    if (err_o < this->orbThrs or this->orbThrs < 0.0) conv_o = true;
    if (err_p < this->propThrs or this->propThrs < 0.0) conv_p = true;
    return (conv_o and conv_p);
}

bool SCF::needLocalization(int nIter) const {
    bool loc = false;
    if (this->canonical) {
        loc = false;
    } else if (nIter <= 2) {
        loc = true;
    } else if (this->rotation == 0) {
        loc = false;
    } else if (nIter%this->rotation == 0) {
        loc = true;
    }
    return loc;
}

bool SCF::needDiagonalization(int nIter) const {
    bool diag = false;
    if (not this->canonical) {
        diag = false;
    } else if (nIter <= 2) {
        diag = true;
    } else if (this->rotation == 0) {
        diag = false;
    } else if (nIter%this->rotation == 0) {
        diag = true;
    }
    return diag;
}

void SCF::printUpdate(const std::string &name, double P, double dP) const {
    int oldPrec = Printer::setPrecision(15);
    double p = 1.0;
    if (std::abs(P) > mrcpp::MachineZero) {
        p = P;
    }
    bool done = (std::abs(dP/p) < this->propThrs) or this->propThrs < 0.0;
    printout(0, name);
    printout(0, std::setw(24) << P);
    Printer::setPrecision(5);
    printout(0, std::setw(16) << dP);
    println(0, std::setw(5) << done);
    Printer::setPrecision(oldPrec);
}

double SCF::getUpdate(const std::vector<double> &vec, int i, bool absPrec) const {
    if (i < 1 or i > vec.size()) MSG_ERROR("Invalid argument");
    double E_i = vec[i-1];
    double E_im1 = 0.0;
    if (i > 1) {
        E_im1 = vec[i-2];
    }
    double E_diff = E_i - E_im1;
    if (not absPrec and std::abs(E_i) > mrcpp::MachineZero) {
        E_diff *= 1.0/E_i;
    }
    return E_diff;
}

void SCF::printOrbitals(const DoubleVector &epsilon, const OrbitalVector &Phi, int flag) const {
    Printer::printHeader(0, "Orbitals");
    if (flag == 0) println(0, " Orb    F(i,i)        Error         nNodes  Spin  Occ  Done ");
    if (flag == 1) println(0, " Orb    Norm          Error         nNodes  Spin  Occ  Done ");
    Printer::printSeparator(0, '-');
    int oldprec = Printer::setPrecision(5);
    bool tot_conv = true;
    for (int i = 0; i < Phi.size(); i++) {
        bool converged = (Phi[i].error() < this->orbThrs) ? true : false;
        printout(0, std::setw(3) << i);
        printout(0, " " << std::setw(13) << epsilon(i));
        printout(0, " " << std::setw(13) << Phi[i].error());
        printout(0, " " << std::setw(10) << Phi[i].getNNodes());
        printout(0, std::setw(5) << Phi[i].printSpin());
        printout(0, std::setw(5) << Phi[i].occ());
        printout(0, std::setw(5) << converged << std::endl);
        if (not converged) tot_conv = false;
    }

    DoubleVector errors = orbital::get_errors(Phi);
    double tot_error = std::sqrt(errors.dot(errors));

    Printer::printSeparator(0, '-');
    printout(0, " Total error:                    ");
    printout(0, std::setw(19) << tot_error << "  ");
    printout(0, std::setw(3) << tot_conv << std::endl);
    Printer::printSeparator(0, '=', 2);
    Printer::setPrecision(oldprec);
}

void SCF::printConvergence(bool converged) const {
    int iter = this->orbError.size();
    int oldPrec = Printer::getPrecision();
    Printer::printHeader(0, "Convergence rate");
    println(0,"Iter    OrbError       Property                   Update  ");
    Printer::printSeparator(0, '-');
    for (int i = 0; i < iter; i++) {
        double prop_i = this->property[i];
        double propDiff = getUpdate(this->property, i+1, true);
        printout(0, std::setw(3) << i+1);
        Printer::setPrecision(5);
        printout(0, std::setw(15) << this->orbError[i]);
        Printer::setPrecision(15);
        printout(0, std::setw(26) << prop_i);
        Printer::setPrecision(5);
        printout(0, std::setw(15) << propDiff);
        printout(0, std::endl);
    }
    Printer::setPrecision(oldPrec);
    Printer::printSeparator(0, '-');
    if (converged) {
        println(0,"                      SCF converged!!!                      ");
    } else {
        println(0,"                   SCF did NOT converge!!!                  ");
    }
    Printer::printSeparator(0, '=', 2);
}

void SCF::printCycle(int nIter) const {
    printout(0, std::endl << std::endl);
    printout(0, "#######################");
    printout(0, " SCF cycle " << std::setw(2) << nIter << " ");
    printout(0, "#######################");
    printout(0, std::endl << std::endl << std::endl);
}

void SCF::printTimer(double t) const {
    int oldPrec = Printer::setPrecision(5);
    printout(0, std::endl << std::endl);
    printout(0, "################");
    printout(0, " Wall time: " << t << " sec ");
    printout(0, "################");
    printout(0, std::endl << std::endl << std::endl);
    Printer::setPrecision(oldPrec);
}

} //namespace mrchem
