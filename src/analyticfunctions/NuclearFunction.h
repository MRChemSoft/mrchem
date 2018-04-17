#pragma once

#include <vector>

#include "MRCPP/MWFunctions"

#include "Nucleus.h"

namespace mrchem {

class NuclearFunction final : public mrcpp::RepresentableFunction<3> {
public:
    NuclearFunction() { }
    ~NuclearFunction() { }

    double evalf(const double *r) const;

    void push_back(const Nucleus &nuc, double S);

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

protected:
    Nuclei nuclei;
    std::vector<double> smooth;
};

} //namespace mrchem
