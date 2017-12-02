#pragma once

#include "QMPotential.h"
#include "Nucleus.h"
#include "NuclearFunction.h"
#include "MWProjector.h"

class NuclearPotential : public QMPotential {
public:
    NuclearPotential(double Z, const double *R = 0, double S = 1.0e-7) {
        func.push_back(Z, R, S);
    }
    NuclearPotential(const Nuclei &nucs, double prec)
        : nuclei(nucs), func(nucs, prec) {
    }
    virtual ~NuclearPotential() { }

    virtual void setup(double prec) {
        if (IS_EQUAL(prec, this->apply_prec)) return;

        setApplyPrec(prec);
        if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
        if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

        MWProjector<3> project(this->apply_prec, this->max_scale);

        Timer timer;
        this->allocReal();
        project(this->real(), this->func);
        timer.stop();

        int n = this->getNNodes();
        double t = timer.getWallTime();
        TelePrompter::printTree(0, "Nuclear potential", n, t);
    }
    virtual void clear() {
        clearReal(true);
        clearImag(true);
        clearApplyPrec();
    }

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

protected:
    Nuclei nuclei;
    NuclearFunction func;
};

