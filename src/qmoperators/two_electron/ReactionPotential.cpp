#include <cmath>
#include "MRCPP/MWOperators"
#include "MRCPP/Gaussians"

#include "ReactionPotential.h"
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/Density.h"
#include "chemistry/chemistry_utils.h"

using mrcpp::FunctionTree;
using mrcpp::PoissonOperator;
using mrcpp::ABGVOperator;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

ReactionPotential::ReactionPotential(mrcpp::PoissonOperator *P, mrcpp::DerivativeOperator<3> *D, Cavity *C, const Nuclei &nucs, OrbitalVector *Phi)
    : QMPotential(1, false)
    , cavity(C)
    , nuclei(nucs)
    , orbitals(Phi)
    , poisson(P)
    , derivative(D)
    , rho_tot(false)
    , rho_el(false)
    , rho_nuc(false)
    , V_n_func(false) {
    //, cavity_func(false)
    //, inv_eps_func(false)
    //, gamma_func(false)
    //, rho_eff_func(false)
    this->e_Energy = 0.0;
    this->nuc_Energy = 0.0;
    this->tot_Energy = 0.0;
}

//  ~ReactionPotential();

//kun en om gangen
QMFunction ReactionPotential::setup_eps(bool is_eps) {
  QMFunction cavity_func;
  cavity_func.alloc(NUMBER::Real);

  cavity->eval_epsilon(is_eps, cavity->islinear());
  qmfunction::project(cavity_func, *cavity, NUMBER::Real, this->apply_prec);

  return cavity_func;
}



QMFunction ReactionPotential::calc_rho_eff(QMFunction inv_eps_func) {
  QMFunction rho_eff_func;
  Density rho_n = chemistry::compute_nuclear_density(this->apply_prec, this->nuclei, 1000);
  Density rho_e(false);
  Density rho(false);
  density::compute(this->apply_prec, rho_e, *orbitals, DENSITY::Total);
  qmfunction::add(rho, -1.0, rho_e, 1.0, rho_n, -1.0);

  this->rho_tot = rho;
  this->rho_el  = rho_e;
  this->rho_nuc = rho_n;
  this->rho_el.rescale(-1.0);
  qmfunction::multiply(rho_eff_func, rho_tot, inv_eps_func, this->apply_prec);
  return rho_eff_func;
}


QMFunction ReactionPotential::calc_gamma(QMFunction inv_eps_func, mrcpp::FunctionTreeVector<3> d_cavity) {
    QMFunction gamma_func;
    gamma_func.alloc(NUMBER::Real);

  auto d_V_n = mrcpp::gradient(*derivative, V_n_func.real());

  if(cavity->islinear()){
    QMFunction temp_func;
    temp_func.alloc(NUMBER::Real);
    mrcpp::dot(this->apply_prec, temp_func.real(), d_V_n, d_cavity);
    qmfunction::multiply(gamma_func, temp_func, inv_eps_func, this->apply_prec);

  }else{
    mrcpp::dot(this->apply_prec, gamma_func.real(), d_V_n, d_cavity);
  }
  gamma_func.rescale(1.0/(4.0*MATHCONST::pi));

  return gamma_func;
}


void ReactionPotential::setup(double prec) {
  setApplyPrec(prec);

  QMFunction cavity_func;
  QMFunction inv_eps_func;
  QMFunction rho_eff_func;
  mrcpp::FunctionTreeVector<3> d_cavity;

  cavity_func  = setup_eps(false);
  inv_eps_func = setup_eps(true);
  cavity_func.rescale(cavity->dcoeff);
  rho_eff_func = calc_rho_eff(inv_eps_func);
  d_cavity = mrcpp::gradient(*derivative, cavity_func.real());

  if(not V_n_func.hasReal()){
    V_n_func.alloc(NUMBER::Real);
    mrcpp::apply(prec, V_n_func.real(), *poisson, rho_eff_func.real());
  }

  double error = 1;
  int i = 1;
  while(error >= this->apply_prec) {
    QMFunction gamma_func;
    QMFunction temp_func;
    QMFunction V_np1_func;
    QMFunction diff_func;

    gamma_func = calc_gamma(inv_eps_func, d_cavity);
    V_np1_func.alloc(NUMBER::Real);

    qmfunction::add(temp_func, 1.0, rho_eff_func, 1.0, gamma_func, -1.0);
    mrcpp::apply(prec, V_np1_func.real(), *poisson, temp_func.real());

    qmfunction::add(diff_func, 1.0, V_n_func, -1.0, V_np1_func, -1.0);
    error = diff_func.norm();
    V_n_func = V_np1_func;
    std::cout << error << ' ' << i << std::endl;
    i++;
  }
  QMFunction V_0_func;

  V_0_func.alloc(NUMBER::Real);
  mrcpp::apply(prec, V_0_func.real(), *poisson, rho_tot.real());

  qmfunction::add(*this, 1.0, V_n_func, -1.0, V_0_func, -1.0);

  cavity_func.free(NUMBER::Real);
  inv_eps_func.free(NUMBER::Real);
  rho_eff_func.free(NUMBER::Real);
  V_0_func.free(NUMBER::Real);

}

double &ReactionPotential::get_tot_Energy(){
  QMFunction temp_prod_func;
  qmfunction::multiply(temp_prod_func, rho_tot, *this, this->apply_prec);
  tot_Energy = temp_prod_func.integrate().real();
  return tot_Energy;
}


double &ReactionPotential::get_e_Energy(){
  QMFunction temp_prod_func;
  qmfunction::multiply(temp_prod_func, rho_el, *this, this->apply_prec);
  e_Energy = temp_prod_func.integrate().real();
  return e_Energy;
}


double &ReactionPotential::get_nuc_Energy(){
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

} //namespace mrchem
