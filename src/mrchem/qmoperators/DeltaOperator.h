#pragma once

#include "QMPotential.h"
#include "GaussFunc.h"

class DeltaOperator : public QMPotential {
public:
    DeltaOperator(const double *o, double expo)
            : QMPotential(1) {
        setPosition(this->r_O, o);
        setFunction(expo, this->r_O);
    };
    virtual ~DeltaOperator() { }

    virtual void setup(double prec);
    virtual void clear();

protected:
    double r_O[3];
    GaussFunc<3> func;

    void setPosition(double *out, const double *inp);
    void setFunction(double beta, const double *pos);
};

