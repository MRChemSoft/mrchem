#include "ReactionPotential.h"
#include "MRCPP/MWOperators"
#include "MRCPP/Plotter"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using Cavity_p = std::shared_ptr<mrchem::Cavity>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

ReactionPotential::ReactionPotential(PoissonOperator_p P, DerivativeOperator_p D, int hist)
        : QMPotential(1, false)
        , history(hist)
        , poisson(P)
        , derivative(D) {}

void ReactionPotential::accelerateConvergence(QMFunction &diff_func, QMFunction &temp, KAIN &kain) {
    OrbitalVector phi_n(0);
    OrbitalVector dPhi_n(0);
    phi_n.push_back(Orbital(SPIN::Paired));
    dPhi_n.push_back(Orbital(SPIN::Paired));

    phi_n[0].QMFunction::operator=(temp);
    dPhi_n[0].QMFunction::operator=(diff_func);

    kain.accelerate(this->apply_prec, phi_n, dPhi_n);

    temp.QMFunction::operator=(phi_n[0]);
    diff_func.QMFunction::operator=(dPhi_n[0]);

    phi_n.clear();
    dPhi_n.clear();
}

void ReactionPotential::setup(double prec) {
    setApplyPrec(prec);
    QMFunction &temp = *this;

    // Solve the poisson equation
    QMFunctionVector Terms = this->helper->makeTerms(this->apply_prec);
    QMFunction V_nm1 = Terms[0];
    QMFunction gamma = Terms[1];
    QMFunction rho_eff = Terms[2];
    QMFunction V_vac = Terms[3];

    if (this->variational) {
        QMFunction poisson_func;
        QMFunction V_n;
        QMFunction dV_n;
        poisson_func.alloc(NUMBER::Real);
        V_n.alloc(NUMBER::Real);
        dV_n.alloc(NUMBER::Real);

        qmfunction::add(poisson_func, 1.0, gamma, 1.0, rho_eff, -1.0);
        mrcpp::apply(this->apply_prec, V_n.real(), *poisson, poisson_func.real());
        qmfunction::add(dV_n, 1.0, V_n, -1.0, V_nm1, -1.0);
        this->helper->updateDifferencePotential(dV_n);
        auto error = dV_n.norm();
        temp = V_n;
        println(0, "error:");
        println(0, error);
    } else {
        KAIN kain(this->history);
        QMFunction V_n;
        QMFunction dV_n;
        double error = 10;
        for (int iter = 1; error >= this->apply_prec; iter++) {
            if (iter > 100) break;
            QMFunction poisson_func;
            QMFunction V_tot;
            this->helper->resetQMFunction(V_n);
            this->helper->resetQMFunction(dV_n);

            // solve the poisson equation
            qmfunction::add(poisson_func, 1.0, gamma, 1.0, rho_eff, -1.0);
            mrcpp::apply(this->apply_prec, V_n.real(), *poisson, poisson_func.real());
            qmfunction::add(dV_n, 1.0, V_n, -1.0, V_nm1, -1.0);

            // use a convergence accelerator
            if (iter > 1 and this->history > 0) accelerateConvergence(dV_n, V_nm1, kain);
            V_n.free(NUMBER::Real);
            qmfunction::add(V_n, 1.0, V_nm1, 1.0, dV_n, -1.0);
            error = dV_n.norm();

            // set up for next iteration
            qmfunction::add(V_tot, 1.0, V_n, 1.0, V_vac, -1.0);
            gamma.free(NUMBER::Real);
            gamma = this->helper->updateGamma(V_tot, this->apply_prec);

            V_nm1.free(NUMBER::Real);
            qmfunction::deep_copy(V_nm1, V_n);

            println(0, "error:");
            println(0, error);
            println(0, "Microiteration:") println(0, iter);
            this->helper->resetQMFunction(poisson_func);
            this->helper->resetQMFunction(V_tot);
        }
        this->helper->updatePotential(V_nm1);
        this->helper->updateDifferencePotential(dV_n);
        temp = V_n;
    }
    V_nm1.free(NUMBER::Real);
    gamma.free(NUMBER::Real);
    rho_eff.free(NUMBER::Real);
    V_vac.free(NUMBER::Real);
    Terms.clear(); 

}

void ReactionPotential::updatePotential(QMFunction new_potential) {
    QMFunction &temp = *this;
    temp.free(NUMBER::Real);
    temp = new_potential;
}

void ReactionPotential::clear() {
    QMFunction::free(NUMBER::Real);
    clearApplyPrec();
    this->helper->clear();
}

} // namespace mrchem
