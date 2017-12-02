#pragma once

#include "QMPotential.h"
#include "NuclearFunction.h"

class HSFPotential : public QMPotential {
public:
    HSFPotential(const Nuclei &nucs, const double *o) : nuclei(nucs) {
        setPosition(this->r_O, o);
    }
    virtual ~HSFPotential() { }

    virtual void setup(double prec) {
        if (IS_EQUAL(prec, this->apply_prec)) return;

        setApplyPrec(prec);
        if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
        if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

        const double *r_O = this->r_O;
        Nuclei &nuclei = this->nuclei;

        auto r_m1 = [] (const double *R, const double *r) -> double {
            NuclearFunction nuc_func;
            nuc_func.push_back(1.0, R, 1.0e-7);
            return nuc_func.evalf(r);
        };

        auto r_m3 = [] (const double *R, const double *r) -> double {
            NuclearFunction nuc_func;
            nuc_func.push_back(1.0, R, 1.0e-7);

            double r_m1 = nuc_func.evalf(r);
            return r_m1*r_m1*r_m1;
        };

        auto rO_dot_rA = [] (const double *r_i, const double *r_O, const double *r_A) -> double {
            double x_O = r_i[0] - r_O[0];
            double y_O = r_i[1] - r_O[1];
            double z_O = r_i[2] - r_O[2];
            double x_A = r_i[0] - r_A[0];
            double y_A = r_i[1] - r_A[1];
            double z_A = r_i[2] - r_A[2];
            return (x_O*x_A + y_O*y_A + z_O*z_A);
        };

        auto f = [r_O, nuclei, r_m1, r_m3, rO_dot_rA] (const double *r) -> double {
            double O_m1 = r_m1(r_O, r);
            double result = 0.0;
            for (int A = 0; A < nuclei.size(); A++) {
                const double *R_A = nuclei[A].getCoord();
                double Z_A = nuclei[A].getCharge();
                double A_m3 = r_m3(R_A, r);
                double O_dot_A = rO_dot_rA(r, r_O, R_A);
                result += Z_A*O_m1*A_m3*O_dot_A;
            }
            return result;
        };

        MWProjector<3> project(this->apply_prec, this->max_scale);

        Timer timer;
        this->allocReal();
        project(this->real(), f);
        timer.stop();

        int n = this->getNNodes();
        double t = timer.getWallTime();
        TelePrompter::printTree(0, "Cubic potential", n, t);
    }

    virtual void clear() {
        clearReal(true);
        clearImag(true);
        clearApplyPrec();
    }

protected:
    double r_O[3];
    Nuclei nuclei;

    void setPosition(double *out, const double *inp) {
        if (inp != 0) {
            out[0] = inp[0];
            out[1] = inp[1];
            out[2] = inp[2];
        } else {
            out[0] = 0.0;
            out[1] = 0.0;
            out[2] = 0.0;
        }
    }
};

