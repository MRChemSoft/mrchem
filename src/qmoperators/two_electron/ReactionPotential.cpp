#include "ReactionPotential.h"
#include "MRCPP/MWOperators"
#include "MRCPP/Plotter"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"
#include <string>

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

    if (this->run_once) {
        this->run_once = false;
        return;
    }
    QMFunction &temp = *this;

    // Solve the poisson equation
    QMFunctionVector Terms = this->helper->makeTerms(this->derivative, this->poisson, this->apply_prec);
    QMFunction V_nm1 = Terms[0];
    QMFunction gamma = Terms[1];
    QMFunction rho_eff = Terms[2];
    QMFunction V_vac = Terms[3];
    QMFunction V_n;
    QMFunction dV_n;
    QMFunction poisson_func;
    if (this->variational) {
        V_n.alloc(NUMBER::Real);
        qmfunction::add(poisson_func, 1.0, gamma, 1.0, rho_eff, -1.0);
        mrcpp::apply(this->apply_prec, V_n.real(), *poisson, poisson_func.real());
        qmfunction::add(dV_n, 1.0, V_n, -1.0, V_nm1, -1.0);
        auto error = dV_n.norm();

        print_utils::text(0, "error:           ", print_utils::dbl_to_str(error, 5, true));
    } else {
        print_utils::headline(0, "Calculating Reaction Potential");
        QMFunction V_tot;

        KAIN kain(this->history);
        double error = 10;
        double converge_prec = this->mo_residual;
        for (int iter = 1; error >= converge_prec && iter <= 100; iter++) {
            if (iter > 1) {
                dV_n.free(NUMBER::Real);
                poisson_func.free(NUMBER::Real);
                V_tot.free(NUMBER::Real);
                V_nm1.free(NUMBER::Real);
                qmfunction::deep_copy(V_nm1, V_n);
            }
            this->helper->resetQMFunction(V_n);

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
            this->helper->updateGamma(gamma, derivative, V_tot, this->apply_prec);

            print_utils::text(0, "error:           ", print_utils::dbl_to_str(error, 5, true));
            print_utils::text(0, "Microiteration:  ", std::to_string(iter));
            // if ((error <= converge_prec) && (iter < 20) && (converge_prec > this->apply_prec)) converge_prec /= 10.0;
        }
        println(0, " Converged Reaction Potential!");
        V_tot.free(NUMBER::Real);
        kain.clear();
    }

    this->helper->updateDifferencePotential(dV_n);
    qmfunction::deep_copy(temp, V_nm1);
    V_nm1.free(NUMBER::Real);
    V_n.free(NUMBER::Real);
    dV_n.free(NUMBER::Real);
    gamma.free(NUMBER::Real);
    rho_eff.free(NUMBER::Real);
    V_vac.free(NUMBER::Real);
    poisson_func.free(NUMBER::Real);
    Terms.clear();
}

void ReactionPotential::clear() {
    //    QMFunction::free(NUMBER::Real);
    clearApplyPrec();
    this->helper->clear();
}

} // namespace mrchem
