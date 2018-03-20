#pragma once

#include "MRCPP/MWFunctions"

namespace mrchem {

class RadialFunction final : public mrcpp::RepresentableFunction<1> {
public:
    RadialFunction(int n, int l, double Z);
    ~RadialFunction() { }

    double evalf(const double *r) const;

protected:
    const int N;
    const int L;
    double c_0;
    double c_1;

    double calcConstant(double Z) const;
    double evalfPoly(double r) const;
};


class AngularFunction final : public mrcpp::RepresentableFunction<3> {
public:
    AngularFunction(int l, int m);
    ~AngularFunction() { }

    double evalf(const double *r) const;

protected:
    const int L;
    const int M;
    double c_0;

    double calcConstant() const;
    double evalfPoly(const double *q) const;
};

class HydrogenFunction final : public mrcpp::RepresentableFunction<3> {
public:
    HydrogenFunction(int n, int l, int m, double Z = 1.0, const double *o = 0);
    ~HydrogenFunction() { }

    double evalf(const double *p) const;

protected:
    double origin[3];
    RadialFunction R;
    AngularFunction Y;
};

} //namespace mrchem
