#pragma once

#include "PositionOperator.h"
#include "Nucleus.h"

namespace mrchem {

class H_E_dip final : public PositionOperator {
public:
    H_E_dip(const double *o = 0) : PositionOperator(o) { 
        RankOneTensorOperator<3> &h = (*this);
        h[0] = (*this)[0];
        h[1] = (*this)[1];
        h[2] = (*this)[2];
        h[0] = -1.0 * h[0];
        h[1] = -1.0 * h[1];
        h[2] = -1.0 * h[2];
    }
    ~H_E_dip() { }

    ComplexVector trace(const Nuclei &nucs) {
        ComplexVector result = ComplexVector::Zero(3);
        for (int k = 0; k < nucs.size(); k++) {
            result += trace(nucs[k]);
        }
        return result;
    }
    ComplexVector trace(const Nucleus &nuc) {
        ComplexVector result = ComplexVector::Zero(3);
        double Z = nuc.getCharge();
        const double *R = nuc.getCoord();
        result(0) = Z*this->r_x.real().evalf(R);
        result(1) = Z*this->r_y.real().evalf(R);
        result(2) = Z*this->r_z.real().evalf(R);
        return result;
    }

    using RankOneTensorOperator<3>::trace;
};

} //namespace mrchem
