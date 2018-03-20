#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"
#include "NuclearFunction.h"

namespace mrchem {

class QMNucPot final : public QMPotential {
public:
    QMNucPot(const Nuclei &nucs, double prec);
    ~QMNucPot() { }

    void setup(double prec);
    void clear();

    Nuclei &getNuclei() { return this->func.getNuclei(); }
    const Nuclei &getNuclei() const { return this->func.getNuclei(); }

protected:
    NuclearFunction func;
};

class NuclearPotential final : public RankZeroTensorOperator {
public:
    NuclearPotential(const Nuclei &nucs, double prec)
            : r_m1(nucs, prec) {
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
    }
    ~NuclearPotential() { }

protected:
    QMNucPot r_m1;
};

} //namespace mrchem
