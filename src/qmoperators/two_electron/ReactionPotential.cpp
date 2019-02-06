#include "MRCPP/Gaussians"
#include "MRCPP/MWOperators"
#include <cmath>

#include "ReactionPotential.h"
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::ABGVOperator;
using mrcpp::FunctionTree;
using mrcpp::PoissonOperator;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

ReactionPotential::ReactionPotential(mrcpp::PoissonOperator *P,
                                     mrcpp::DerivativeOperator<3> *D,
                                     Cavity *C,
                                     const Nuclei &nucs,
                                     OrbitalVector *Phi,
                                     bool testing)
        : QMPotential(1, false)
        , cavity(C)
        , nuclei(nucs)
        , orbitals(Phi)
        , poisson(P)
        , derivative(D)
        , rho_tot(false)
        , rho_el(false)
        , rho_nuc(false) {
    this->e_Energy = 0.0;
    this->nuc_Energy = 0.0;
    this->tot_Energy = 0.0;
    this->testing = testing;
}

void ReactionPotential::calc_eps(bool is_inv, QMFunction &cavity_func) {
    cavity->eval_epsilon(is_inv, cavity->islinear());
    qmfunction::project(cavity_func, *cavity, NUMBER::Real, this->apply_prec);
}

void ReactionPotential::calc_rho_eff(QMFunction const &inv_eps_func,
                                     QMFunction &rho_eff_func,
                                     QMFunction &cavity_func) {
    rho_nuc = chemistry::compute_nuclear_density(this->apply_prec, this->nuclei, 1000);
    density::compute(this->apply_prec, rho_el, *orbitals, DENSITY::Total);
    rho_el.rescale(-1.0);
    qmfunction::add(rho_tot, 1.0, rho_el, 1.0, rho_nuc, -1.0);

    QMFunction ones;
    ones.alloc(NUMBER::Real);

    cavity->change_use(true);
    calc_eps(false, ones);
    cavity->change_use(false);

    qmfunction::multiply(rho_eff_func, rho_tot, ones, this->apply_prec);

    ones.free(NUMBER::Real);
}

void ReactionPotential::calc_gamma(QMFunction const &inv_eps_func,
                                   QMFunction &gamma_func,
                                   QMFunction &V_0_func,
                                   QMFunction &temp,
                                   mrcpp::FunctionTreeVector<3> &d_cavity) {
    QMFunction temp_func1;
    QMFunction temp_func2;
    temp_func1.alloc(NUMBER::Real);
    temp_func2.alloc(NUMBER::Real);

    qmfunction::add(temp_func1, 1.0, temp, 1.0, V_0_func, -1.0);

    auto d_V = mrcpp::gradient(*derivative, temp_func1.real());

    mrcpp::dot(this->apply_prec, temp_func2.real(), d_V, d_cavity);
    qmfunction::multiply(gamma_func, temp_func2, inv_eps_func, this->apply_prec);

    gamma_func.rescale(1.0 / (4.0 * MATHCONST::pi));
    temp_func1.free(NUMBER::Real);
    temp_func2.free(NUMBER::Real);
    d_V.clear();
}

void ReactionPotential::grad_G(QMFunction &gamma_func,
                               QMFunction &cavity_func,
                               QMFunction &rho_tot,
                               QMFunction &grad_G_func) {
    QMFunction temp_func;
    qmfunction::multiply(temp_func, gamma_func, cavity_func, this->apply_prec);
    qmfunction::add(grad_G_func, 1.0, temp_func, -1.0, rho_tot, -1.0);
}

void ReactionPotential::setup(double prec) {
    setApplyPrec(prec);

    QMFunction V_func;
    QMFunction V_0_func;
    QMFunction cavity_func;
    QMFunction inv_eps_func;
    QMFunction rho_eff_func;
    QMFunction gamma_func;
    QMFunction &temp = *this;
    mrcpp::FunctionTreeVector<3> d_cavity;

    calc_eps(false, cavity_func);
    calc_eps(true, inv_eps_func);
    calc_rho_eff(inv_eps_func, rho_eff_func, cavity_func);
    d_cavity = mrcpp::gradient(*derivative, cavity_func.real());

    V_0_func.alloc(NUMBER::Real);
    mrcpp::apply(prec, V_0_func.real(), *poisson, rho_tot.real());

    if (not(temp).hasReal()) {
        V_func.alloc(NUMBER::Real);
        mrcpp::apply(prec, V_func.real(), *poisson, rho_eff_func.real());
        qmfunction::add(temp, 1.0, V_func, -1.0, V_0_func, -1.0);
    }

    auto error = 1.00;
    int iter = 0;
    while (error >= this->apply_prec) {
        gamma_func.free(NUMBER::Real);
        QMFunction temp_func;
        QMFunction V_np1_func;
        QMFunction diff_func;

        calc_gamma(inv_eps_func, gamma_func, V_0_func, temp, d_cavity);
        V_np1_func.alloc(NUMBER::Real);

        qmfunction::add(temp_func, 1.0, rho_eff_func, 1.0, gamma_func, -1.0);
        mrcpp::apply(prec, V_np1_func.real(), *poisson, temp_func.real());

        qmfunction::add(diff_func, 1.0, temp, -1.0, V_np1_func, -1.0);
        error = diff_func.norm();

        temp = V_np1_func;
        iter++;
        std::cout << "iter.:\t" << iter << "\n"
                  << "error:\t" << error << std::endl;
        if (this->testing and iter > 3) break;
        if (error >= 100000.00) break;
    }

    /*
        V_0_func.alloc(NUMBER::Real);
        mrcpp::apply(prec, V_0_func.real(), *poisson, rho_tot.real());

        qmfunction::add(*this, 1.0, V_func, -1.0, V_0_func, -1.0);*/

    d_cavity.clear();
    gamma_func.free(NUMBER::Real);
    cavity_func.free(NUMBER::Real);
    inv_eps_func.free(NUMBER::Real);
    rho_eff_func.free(NUMBER::Real);
    V_0_func.free(NUMBER::Real);
    V_func.free(NUMBER::Real);
}

double &ReactionPotential::get_tot_Energy() {
    QMFunction temp_prod_func;
    qmfunction::multiply(temp_prod_func, rho_tot, *this, this->apply_prec);
    tot_Energy = temp_prod_func.integrate().real();
    return tot_Energy;
}

double &ReactionPotential::get_e_Energy() {
    QMFunction temp_prod_func;
    qmfunction::multiply(temp_prod_func, rho_el, *this, this->apply_prec);
    e_Energy = temp_prod_func.integrate().real();
    return e_Energy;
}

double &ReactionPotential::get_nuc_Energy() {
    QMFunction temp_prod_func;
    qmfunction::multiply(temp_prod_func, rho_nuc, *this, this->apply_prec);
    nuc_Energy = temp_prod_func.integrate().real();
    return nuc_Energy;
}

void ReactionPotential::clear() {
    clearApplyPrec();
    rho_tot.free(NUMBER::Real);
    rho_el.free(NUMBER::Real);
    rho_nuc.free(NUMBER::Real);
    QMFunction::free(NUMBER::Total);
}

} // namespace mrchem
