#pragma once

#include <cmath>

#include "MRCPP/MWFunctions"

namespace mrchem {

class RadialFunction final : public mrcpp::RepresentableFunction<1> {
public:
    RadialFunction(int n, int l, double Z)
            : N(n), L(l) {
        this->c_0 = calcConstant(Z);
        this->c_1 = 2.0*Z/(1.0*this->N);
    }
    ~RadialFunction() { }

    double evalf(const double *r) const {
        double rho = this->c_1*r[0];
        return this->c_0*this->evalfPoly(rho)*exp(-rho/2.0);
    }

protected:
    const int N;
    const int L;
    double c_0;
    double c_1;

    double calcConstant(double Z) const {
        double c = 0.0;
               if (N == 1 and L == 0) { c = 2.0;
        } else if (N == 2 and L == 0) { c = 1.0/(   2.0*std::sqrt(  2.0));
        } else if (N == 3 and L == 0) { c = 1.0/(   9.0*std::sqrt(  3.0));
        } else if (N == 4 and L == 0) { c = 1.0/(  48.0*std::sqrt(  4.0));
        } else if (N == 5 and L == 0) { c = 1.0/( 300.0*std::sqrt(  5.0));
        } else if (N == 6 and L == 0) { c = 1.0/(2160.0*std::sqrt(  6.0));
        } else if (N == 2 and L == 1) { c = 1.0/(   2.0*std::sqrt(  6.0));
        } else if (N == 3 and L == 1) { c = 1.0/(   9.0*std::sqrt(  6.0));
        } else if (N == 4 and L == 1) { c = 1.0/(  32.0*std::sqrt( 15.0));
        } else if (N == 5 and L == 1) { c = 1.0/( 150.0*std::sqrt( 30.0));
        } else if (N == 6 and L == 1) { c = 1.0/( 432.0*std::sqrt(210.0));
        } else if (N == 3 and L == 2) { c = 1.0/(   9.0*std::sqrt( 30.0));
        } else if (N == 4 and L == 2) { c = 1.0/(  96.0*std::sqrt(  5.0));
        } else if (N == 5 and L == 2) { c = 1.0/( 150.0*std::sqrt( 70.0));
        } else if (N == 6 and L == 2) { c = 1.0/( 864.0*std::sqrt(105.0));
        } else if (N == 4 and L == 3) { c = 1.0/(  96.0*std::sqrt( 35.0));
        } else if (N == 5 and L == 3) { c = 1.0/( 300.0*std::sqrt( 70.0));
        } else { NOT_IMPLEMENTED_ABORT; }
        return c*pow(Z, 3.0/2.0);
    }

    double evalfPoly(double r) const {
        double value = 0.0;
               if (N == 1 and L == 0) { value =        (  1.0                                                                     );
        } else if (N == 2 and L == 0) { value =        (  2.0 -    1.0*r                                                          );
        } else if (N == 3 and L == 0) { value =        (  6.0 -    6.0*r +    1.0*r*r                                             );
        } else if (N == 4 and L == 0) { value =        ( 24.0 -   36.0*r +   12.0*r*r -   1.0*r*r*r                               );
        } else if (N == 5 and L == 0) { value =        (120.0 -  240.0*r +  120.0*r*r -  20.0*r*r*r +  1.0*r*r*r*r                );
        } else if (N == 6 and L == 0) { value =        (720.0 - 1800.0*r + 1200.0*r*r - 300.0*r*r*r + 30.0*r*r*r*r - 1.0*r*r*r*r*r);
        } else if (N == 2 and L == 1) { value =      r*(  1.0                                                                     );
        } else if (N == 3 and L == 1) { value =      r*(  4.0 -    1.0*r                                                          );
        } else if (N == 4 and L == 1) { value =      r*( 20.0 -   10.0*r +    1.0*r*r                                             );
        } else if (N == 5 and L == 1) { value =      r*(120.0 -   90.0*r +   18.0*r*r -   1.0*r*r*r                               );
        } else if (N == 6 and L == 1) { value =      r*(840.0 -  840.0*r +  252.0*r*r -  28.0*r*r*r +  1.0*r*r*r*r                );
        } else if (N == 3 and L == 2) { value =    r*r*(  1.0                                                                     );
        } else if (N == 4 and L == 2) { value =    r*r*(  6.0 -    1.0*r                                                          );
        } else if (N == 5 and L == 2) { value =    r*r*( 42.0 -   14.0*r +    1.0*r*r                                             );
        } else if (N == 6 and L == 2) { value =    r*r*(336.0 -  168.0*r +   24.0*r*r -   1.0*r*r*r                               );
        } else if (N == 4 and L == 3) { value =  r*r*r*(  1.0                                                                     );
        } else if (N == 5 and L == 3) { value =  r*r*r*(  8.0 -    1.0*r                                                          );
        } else { NOT_IMPLEMENTED_ABORT; }
        return value;
    }
};


class AngularFunction final : public mrcpp::RepresentableFunction<3> {
public:
    AngularFunction(int l, int m)
            : L(l), M(m) {
        this->c_0 = calcConstant();
    }
    ~AngularFunction() { }

    double evalf(const double *r) const {
        return this->c_0*this->evalfPoly(r);
    }

protected:
    const int L;
    const int M;
    double c_0;

    double calcConstant() const {
        double c = 0.0;
               if (L == 0 and M == 0) { c = std::sqrt(  1.0/1.0);
        } else if (L == 1 and M == 0) { c = std::sqrt(  3.0/1.0);
        } else if (L == 1 and M == 1) { c = std::sqrt(  3.0/1.0);
        } else if (L == 1 and M == 2) { c = std::sqrt(  3.0/1.0);
        } else if (L == 2 and M == 0) { c = std::sqrt( 60.0/4.0);
        } else if (L == 2 and M == 1) { c = std::sqrt( 60.0/4.0);
        } else if (L == 2 and M == 2) { c = std::sqrt( 60.0/4.0);
        } else if (L == 2 and M == 3) { c = std::sqrt( 15.0/4.0);
        } else if (L == 2 and M == 4) { c = std::sqrt(  5.0/4.0);
        } else if (L == 3 and M == 0) { c = std::sqrt(  7.0/4.0);
        } else if (L == 3 and M == 1) { c = std::sqrt(  7.0/4.0);
        } else if (L == 3 and M == 2) { c = std::sqrt(  7.0/4.0);
        } else if (L == 3 and M == 3) { c = std::sqrt(105.0/4.0);
        } else if (L == 3 and M == 4) { c = std::sqrt(105.0/4.0);
        } else if (L == 3 and M == 5) { c = std::sqrt(105.0/4.0);
        } else if (L == 3 and M == 6) { c = std::sqrt(105.0/4.0);
        } else { NOT_IMPLEMENTED_ABORT; }
        return c/std::sqrt(4.0*mrcpp::pi);
    }

    double evalfPoly(const double *q) const {
        double x = q[0];
        double y = q[1];
        double z = q[2];
        double r = std::sqrt(x*x + y*y + z*z);
        double value = 0.0;
               if (L == 0 and M == 0) { value = 1.0;
        } else if (L == 1 and M == 0) { value = x/r;
        } else if (L == 1 and M == 1) { value = y/r;
        } else if (L == 1 and M == 2) { value = z/r;
        } else if (L == 2 and M == 0) { value = x*y/(r*r);
        } else if (L == 2 and M == 1) { value = x*z/(r*r);
        } else if (L == 2 and M == 2) { value = y*z/(r*r);
        } else if (L == 2 and M == 3) { value = (x*x - y*y)/(r*r);
        } else if (L == 2 and M == 4) { value = (2.0*z*z - x*x - y*y)/(r*r);
        } else if (L == 3 and M == 0) { value = x*(5.0*x*x-3.0*r*r)/(r*r*r);
        } else if (L == 3 and M == 1) { value = y*(5.0*y*y-3.0*r*r)/(r*r*r);
        } else if (L == 3 and M == 2) { value = z*(5.0*z*z-3.0*r*r)/(r*r*r);
        } else if (L == 3 and M == 3) { value = x*(z*z - y*y)/(r*r*r);
        } else if (L == 3 and M == 4) { value = y*(z*z - x*x)/(r*r*r);
        } else if (L == 3 and M == 5) { value = z*(x*x - y*y)/(r*r*r);
        } else if (L == 3 and M == 6) { value = 2.0*x*y*z/(r*r*r);
        } else { NOT_IMPLEMENTED_ABORT; }
        return value;
    }
};

class HydrogenFunction final : public mrcpp::RepresentableFunction<3> {
public:
    HydrogenFunction(int n, int l, int m, double Z = 1.0, const double *o = 0)
            : RepresentableFunction<3>(),
              R(n, l, Z),
              Y(l, m) {
        for (int d = 0; d < 3; d++) this->origin[0] = (o != 0) ? o[0] : 0.0;
    }
    ~HydrogenFunction() { }

    double evalf(const double *p) const {
        double q[3];
        q[0] = p[0]-this->origin[0];
        q[1] = p[1]-this->origin[1];
        q[2] = p[2]-this->origin[2];
        double r = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
        return this->R.evalf(&r)*this->Y.evalf(q);
    }

protected:
    double origin[3];
    RadialFunction R;
    AngularFunction Y;
};

} //namespace mrchem
