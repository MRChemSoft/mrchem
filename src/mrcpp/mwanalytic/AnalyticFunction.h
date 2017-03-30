#ifndef ANALYTICFUNCTION_H
#define ANALYTICFUNCTION_H

#include "RepresentableFunction.h"

template<int D>
class AnalyticFunction : public RepresentableFunction<D> {
public:
    AnalyticFunction(std::function<double (const double *r)> f,
                     const double *a = 0,
                     const double *b = 0)
        : RepresentableFunction<D>(a, b),
          func(f) { }
    virtual ~AnalyticFunction() { }

    virtual double evalf(const double *r) const {
        double val = 0.0;
        if (not this->outOfBounds(r)) val = this->func(r);
        return val;
    }
protected:
    std::function<double (const double *r)> func;
};

#endif // ANALYTICFUNCTION_H
