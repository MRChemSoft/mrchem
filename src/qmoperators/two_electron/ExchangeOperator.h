#pragma once

#include "QMOperator.h"
#include "OrbitalVector.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

class ExchangeOperator : public QMOperator {
public:
    ExchangeOperator(PoissonOperator &P, OrbitalVector &phi, double x_fac)
            : QMOperator(MRA->getMaxScale()),
              x_factor(x_fac),
              poisson(&P),
              orbitals(&phi),
              screen(true) {
        int nOrbs = phi.size();
        this->tot_norms = Eigen::VectorXd::Zero(nOrbs);
        this->part_norms = Eigen::MatrixXd::Zero(nOrbs, nOrbs);
    }
    virtual ~ExchangeOperator() { }

    virtual void rotate(Eigen::MatrixXd &U) = 0;

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    void setExchangeFactor(double x_fac) { this->x_factor = x_fac; }
    double getExchangeFactor() const { return this->x_factor; }

    void setScreen(bool s) { this->screen = s; }
    bool getScreen() const { return this->screen; }

protected:
    double x_factor;            ///< Exchange factor for hybrid XC functionals
    PoissonOperator *poisson;   ///< Pointer to external object
    OrbitalVector *orbitals;    ///< Pointer to external object

    bool screen;                ///< Apply screening in exchange evaluation
    Eigen::VectorXd tot_norms;  ///< Total norms for use in screening
    Eigen::MatrixXd part_norms; ///< Partial norms for use in screening

    double getScaledPrecision(int i, int j) const {
        double scaled_prec = this->apply_prec;
        if (getScreen()) {
            double tNorm = this->tot_norms(i);
            double pNorm = std::max(this->part_norms(i,j), this->part_norms(j,i));
            if (tNorm > 0.0) scaled_prec *= tNorm/pNorm;
        }
        return scaled_prec;
    }
};

