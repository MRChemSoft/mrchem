#ifndef ANALYTICPOTENTIAL_H
#define ANALYTICPOTENTIAL_H

#include "QMPotential.h"
#include "MWProjector.h"
#include "AnalyticFunction.h"

class AnalyticPotential : public QMPotential {
public:
    AnalyticPotential(const std::function<double (const double *r)> *func_r = 0,
                      const std::function<double (const double *r)> *func_i = 0)
            : real_func(0), imag_func(0) {
        if (func_r != 0) setReal(*func_r);
        if (func_i != 0) setImag(*func_i);
    }
    virtual ~AnalyticPotential() {
        if (this->real_func != 0) delete this->real_func;
        if (this->imag_func != 0) delete this->imag_func;
    }

    void setReal(const std::function<double (const double *r)> &func) {
        if (this->real_func != 0) MSG_ERROR("Function not empty");
        this->real_func = new AnalyticFunction<3>(func);
    }
    void setImag(const std::function<double (const double *r)> &func) {
        if (this->imag_func != 0) MSG_ERROR("Function not empty");
        this->imag_func = new AnalyticFunction<3>(func);
    }

    virtual void setup(double prec) {
        if (IS_EQUAL(prec, this->apply_prec)) return;

        setApplyPrec(prec);
        if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
        if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

        MWProjector<3> project(this->apply_prec, this->max_scale);

        Timer timer;
        if (this->real_func != 0) {
            this->allocReal();
            project(this->real(), *this->real_func);
        }
        if (this->imag_func != 0) {
            this->allocImag();
            project(this->imag(), *this->imag_func);
        }
        timer.stop();
        int n = this->getNNodes();
        double t = timer.getWallTime();
        TelePrompter::printTree(0, "Analytic potential", n, t);
    }
    virtual void clear() {
        clearReal(true);
        clearImag(true);
        clearApplyPrec();
    }

protected:
    AnalyticFunction<3> *real_func;
    AnalyticFunction<3> *imag_func;
};

#endif // ANALYTICPOTENTIAL_H
