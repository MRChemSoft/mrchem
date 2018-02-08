#include "NuclearFunction.h"
#include "Nucleus.h"
#include "MathUtils.h"

using namespace std;

NuclearFunction::NuclearFunction()
        : RepresentableFunction<3>(),
          const_fac(-1.0/(3.0*root_pi)) {
}

NuclearFunction::NuclearFunction(const Nuclei &nucs, double prec)
        : RepresentableFunction<3>(),
          const_fac(-1.0/(3.0*root_pi)) {
    int oldprec = Printer::setPrecision(5);
    Printer::printHeader(0, "Setting up nuclear potential");
    println(0, " Nr  Element         Charge        Precision     Smoothing ");
    Printer::printSeparator(0, '-');

    double c = 0.00435*prec;
    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        double Z = nuc.getCharge();
        double z_5 = pow(Z, 5.0);
        double smooth = pow(c/z_5, 1.0/3.0);
        this->push_back(nuc, smooth);

        stringstream symbol;
        symbol << nuc.getElement().getSymbol();
        symbol << "  ";
        printout(0, setw(3) << i+1 << "     ");
        printout(0, symbol.str()[0] << symbol.str()[1]);
        printout(0, setw(22) << Z);
        printout(0, setw(14) << prec);
        printout(0, setw(14) << smooth << endl);
    }
    Printer::printSeparator(0, '=', 2);
    Printer::setPrecision(oldprec);
}

void NuclearFunction::push_back(const Nucleus &nuc, double S) {
    double Z = nuc.getCharge();
    const double *R = nuc.getCoord();
    push_back(Z, R, S);
}

void NuclearFunction::push_back(double Z, const double *R, double S) {
    this->x_coords.push_back(R[0]);
    this->y_coords.push_back(R[1]);
    this->z_coords.push_back(R[2]);
    this->charges.push_back(Z);
    this->smoothParam.push_back(S);
}

double NuclearFunction::evalf(const double *x) const {
    double result = 0.0;
    for (int i = 0; i < this->x_coords.size(); i++) {
        double xyz[3];
        xyz[0] = this->x_coords[i];
        xyz[1] = this->y_coords[i];
        xyz[2] = this->z_coords[i];
        double r1 = MathUtils::calcDistance(3, xyz, x);
        r1 *= 1.0/this->smoothParam[i];
        double r2 = pow(r1, 2.0);
        double partResult = -erf(r1)/r1;
        partResult += this->const_fac*(exp(-r2) + 16.0*exp(-4.0*r2));
        result += this->charges[i]*partResult/this->smoothParam[i];
    }
    return result;
}

bool NuclearFunction::isVisibleAtScale(int scale, int nQuadPts) const {
    double minSmooth = this->smoothParam[0];
    for (int i = 1; i < this->smoothParam.size(); i++) {
        if (this->smoothParam[i] < minSmooth) {
            minSmooth = this->smoothParam[i];
        }
    }
    double stdDeviation = pow(minSmooth, -0.5);
    int visibleScale = int (floor(log2(nQuadPts*5.0*stdDeviation)));
    if (scale < visibleScale) {
        return false;
    }
    return true;
}

bool NuclearFunction::isZeroOnInterval(const double *a, const double *b) const {
    int nNucs = this->x_coords.size();
    int totSplit = 0;
    for (int i = 0; i < nNucs; i++) {
        double x_i = this->x_coords[i];
        double y_i = this->y_coords[i];
        double z_i = this->z_coords[i];
        int split = 1;
        if (a[0] > x_i or b[0] < x_i) split = 0;
        if (a[1] > y_i or b[1] < y_i) split = 0;
        if (a[2] > z_i or b[2] < z_i) split = 0;
        totSplit += split;
    }
    if (totSplit == 0) {
        return true;
    } else {
        return false;
    }
}
