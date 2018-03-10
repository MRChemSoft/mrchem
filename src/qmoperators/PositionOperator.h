#pragma once

#include "AnalyticFunction.h"

#include "QMPotential.h"
#include "RankOneTensorOperator.h"

namespace mrchem {

class QMPosition : public QMPotential {
public:
    QMPosition(int d, const double *o);
    virtual ~QMPosition() { }

protected:
    mrcpp::AnalyticFunction<3> func;

    virtual void setup(double prec);
    virtual void clear();
};

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const double *o = 0)
            : r_x(0, o),
              r_y(1, o),
              r_z(2, o) {
        RankOneTensorOperator &r = *this;
        r[0] = r_x;
        r[1] = r_y;
        r[2] = r_z;
    }
    virtual ~PositionOperator() { }

protected:
    QMPosition r_x;
    QMPosition r_y;
    QMPosition r_z;
};

} //namespace mrchem
