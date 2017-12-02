#pragma once

#include "QMTensorOperator.h"
#include "PositionOperator.h"
#include "AngularMomentumOperator.h"
#include "NuclearPotential.h"
#include "QuadraticPotential.h"
#include "CubicPotential.h"
#include "HSFPotential.h"

class HSFOperatorL : public RankZeroTensorOperator {
public:
    HSFOperatorL(DerivativeOperator<3> &D, const double *o = 0)
            : r_O(o),
              l_O(D, o),
              o_m3(1.0, o) {
        initializeTensorOperator();
    }
    virtual ~HSFOperatorL() { }

protected:
    PositionOperator r_O;
    AngularMomentumOperator l_O;
    CubicPotential o_m3;

    void initializeTensorOperator() {
        RankZeroTensorOperator &l_x = this->l_O[0];
        RankZeroTensorOperator &l_y = this->l_O[1];
        RankZeroTensorOperator &l_z = this->l_O[2];
        RankZeroTensorOperator &o_x = this->r_O[0];
        RankZeroTensorOperator &o_y = this->r_O[1];
        RankZeroTensorOperator &o_z = this->r_O[2];

        RankZeroTensorOperator &h = (*this);
        h = 1.0/(2.0*pi)*o_m3*(l_x*l_x + l_y*l_y + l_z*l_z);
    }
};

class HSFOperatorU : public RankZeroTensorOperator {
public:
    HSFOperatorU(const Nuclei &nucs, const double *o = 0)
            : hsf_pot(nucs, o) {
        initializeTensorOperator();
    }
    virtual ~HSFOperatorU() { }

protected:
    HSFPotential hsf_pot;

    void initializeTensorOperator() {
        RankZeroTensorOperator &h = (*this);
        h = 1.0/(2.0*pi)*hsf_pot;
    }
};

class HSFOperatorV : public RankZeroTensorOperator {
public:
    HSFOperatorV(const double *o = 0) {
        initializeTensorOperator();
    }
    virtual ~HSFOperatorV() { }

protected:
    void initializeTensorOperator() {
        NOT_IMPLEMENTED_ABORT;
    }
};


