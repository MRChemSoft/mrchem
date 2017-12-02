#pragma once

#include <vector>

#include "RepresentableFunction.h"

class Nuclei;
class Nucleus;

class NuclearFunction : public RepresentableFunction<3> {
public:
    NuclearFunction();
    NuclearFunction(const Nuclei &nucs, double prec);
    NuclearFunction &operator=(const NuclearFunction &func) { NOT_IMPLEMENTED_ABORT; }
    virtual ~NuclearFunction() { }

    void push_back(const Nucleus &nuc, double S);
    void push_back(double Z, const double *R, double S);

    double evalf(const double *r) const;

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

protected:
    const double const_fac;
    std::vector<double> smoothParam;
    std::vector<double> charges;
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    std::vector<double> z_coords;
};

