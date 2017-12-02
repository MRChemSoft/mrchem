#pragma once

#include "PositionOperator.h"

class H_E_dip : public PositionOperator {
public:
    H_E_dip(const double *o = 0) : PositionOperator(o) { }
    virtual ~H_E_dip() { }

    Eigen::VectorXd trace(const Nuclei &nucs) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(3);
        for (int k = 0; k < nucs.size(); k++) {
            result += trace(nucs[k]);
        }
        return result;
    }
    Eigen::VectorXd trace(const Nucleus &nuc) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(3);
        double Z = nuc.getCharge();
        const double *R = nuc.getCoord();
        result(0) = -Z*this->r_x.real().evalf(R);
        result(1) = -Z*this->r_y.real().evalf(R);
        result(2) = -Z*this->r_z.real().evalf(R);
        return result;
    }

    using RankOneTensorOperator<3>::trace;
};


