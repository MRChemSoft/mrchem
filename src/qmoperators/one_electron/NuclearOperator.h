#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"
#include "NuclearFunction.h"

namespace mrchem {

class NuclearPotential final : public QMPotential {
public:
    NuclearPotential(const Nuclei &nucs, double prec);
    ~NuclearPotential() { }

    void setup(double prec);
    void clear();

    Nuclei &getNuclei() { return this->func.getNuclei(); }
    const Nuclei &getNuclei() const { return this->func.getNuclei(); }

protected:
    NuclearFunction func;
};

class NuclearOperator final : public RankZeroTensorOperator {
public:
    NuclearOperator(const Nuclei &nucs, double prec)
            : r_m1(nucs, prec) {
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
    }
    ~NuclearOperator() { }

    Nuclei &getNuclei() { return this->r_m1.getNuclei(); }
    const Nuclei &getNuclei() const { return this->r_m1.getNuclei(); }

    double trace(const Nuclei &nucs); 

    using RankZeroTensorOperator::trace;
    
 protected:
    NuclearPotential r_m1;
};

} //namespace mrchem
