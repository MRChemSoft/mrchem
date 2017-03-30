#ifndef POSITIONOPERATOR_H
#define POSITIONOPERATOR_H

#include "QMTensorOperator.h"
#include "AnalyticPotential.h"

class QMPosition : public AnalyticPotential {
public:
    QMPosition(int d, const double *o = 0) {
        double orig = 0.0;
        if (o != 0) orig = o[d];
        auto f = [d, orig] (const double *r) -> double { return r[d] - orig; };
        setReal(f);
    }
    virtual ~QMPosition() { }
};

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const double *o = 0) : r_x(0, o), r_y(1, o), r_z(2, o) {
        initializeTensorOperator();
    }
    virtual ~PositionOperator() { }

protected:
    QMPosition r_x;
    QMPosition r_y;
    QMPosition r_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = r_x;
        h[1] = r_y;
        h[2] = r_z;
    }
};

#endif // POSITIONOPERATOR_H
