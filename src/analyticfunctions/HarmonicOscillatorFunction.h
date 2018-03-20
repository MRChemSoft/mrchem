#pragma once

#include "MRCPP/MWFunctions"

#include "mrchem.h"

namespace mrchem {

class HarmonicOscillator1D final {
public:
    HarmonicOscillator1D(int n, double m, double k, double o)
            : nu(n),
              alpha(std::pow(m*k, 1.0/4.0)),
              origin(o) {
    }
    ~HarmonicOscillator1D() { }

    double operator()(double x) const {
        double ax = this->alpha*(x - this->origin);
        double ap = std::sqrt(this->alpha/MATHCONST::sqrt_pi);
        double N_nu = std::sqrt(N2(this->nu));
        return ap*N_nu*H(this->nu, ax)*exp(-ax*ax/2.0);
    }

protected:
    const double nu;
    const double alpha;
    const double origin;

    double N2(int nu) const {
        if (nu == 0) return 1.0;
        return 1.0/(2.0*nu) * N2(nu-1);
    }
    double H(int nu, double x) const {
        if (nu == 0) return 1.0;
        if (nu == 1) return 2.0*x;
        return 2.0*x*H(nu-1, x) - 2.0*(nu-1.0)*H(nu-2, x);
    }
};

class HarmonicOscillatorFunction final : public mrcpp::RepresentableFunction<3> {
public:
    HarmonicOscillatorFunction(int n[3], double m = 1.0, double *k = 0, const double *o = 0)
            : fx(n[0], m, ((k != 0) ? k[0] : 1.0), ((o != 0) ? o[0] : 0.0)),
              fy(n[1], m, ((k != 0) ? k[1] : 1.0), ((o != 0) ? o[1] : 0.0)),
              fz(n[2], m, ((k != 0) ? k[2] : 1.0), ((o != 0) ? o[2] : 0.0)) {
    }
    ~HarmonicOscillatorFunction() { }

    double evalf(const double *r) const { return fx(r[0])*fy(r[1])*fz(r[2]); }

protected:
    const HarmonicOscillator1D fx;
    const HarmonicOscillator1D fy;
    const HarmonicOscillator1D fz;
};

} //namespace mrchem
