#include "MathUtils.h"

#include "NuclearFunction.h"
#include "Nucleus.h"

namespace mrchem {

void NuclearFunction::push_back(const Nucleus &nuc, double S) {
    this->nuclei.push_back(nuc);
    this->smooth.push_back(S);
}

double NuclearFunction::evalf(const double *r) const {
    double c = -1.0/(3.0*mrcpp::root_pi);
    double result = 0.0;
    for (int i = 0; i < this->nuclei.size(); i++) {
        double S_i = this->smooth[i];
        double Z_i = this->nuclei[i].getCharge();
        const double *R = this->nuclei[i].getCoord();
        double R1 = mrcpp::MathUtils::calcDistance(3, R, r)/S_i;
        double partResult = -erf(R1)/R1 + c*(exp(-R1*R1) + 16.0*exp(-4.0*R1*R1));
        result += Z_i*partResult/S_i;
    }
    return result;
}

bool NuclearFunction::isVisibleAtScale(int scale, int nQuadPts) const {
    double minSmooth = this->smooth[0];
    for (int i = 1; i < this->smooth.size(); i++) {
        if (this->smooth[i] < minSmooth) {
            minSmooth = this->smooth[i];
        }
    }
    double stdDeviation = pow(minSmooth, -0.5);
    int visibleScale = int (floor(log2(nQuadPts*5.0*stdDeviation)));
    if (scale < visibleScale) {
        return false;
    } else {
        return true;
    }
}

bool NuclearFunction::isZeroOnInterval(const double *a, const double *b) const {
    int totSplit = 0;
    for (int i = 0; i < this->nuclei.size(); i++) {
        const double *R = this->nuclei[i].getCoord();
        int split = 1;
        if (a[0] > R[0] or b[0] < R[0]) split = 0;
        if (a[1] > R[1] or b[1] < R[1]) split = 0;
        if (a[2] > R[2] or b[2] < R[2]) split = 0;
        totSplit += split;
    }
    if (totSplit == 0) {
        return true;
    } else {
        return false;
    }
}

} //namespace mrchem
