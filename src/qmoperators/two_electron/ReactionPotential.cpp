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
    , cavity_func(false)
    , inv_eps_func(false)
    , rho_eff_func(false)
    , gamma_func(false)
    , V_n_func(false)
    , rho_el(false)
    , rho_nuc(false) {
    this->e_Energy = 0.0;
    this->nuc_Energy = 0.0;
    this->tot_Energy = 0.0;
}

//  ~ReactionPotential();

void ReactionPotential::setup_eps() {
  cavity->change_radius(4.00);
  cavity->eval_epsilon(false, cavity->islinear());

  qmfunction::project(cavity_func, *cavity, NUMBER::Real, this->apply_prec);

  cavity_func.rescale(cavity->dcoeff);

  cavity->eval_epsilon(true, cavity->islinear());

  qmfunction::project(inv_eps_func, *cavity, NUMBER::Real, this->apply_prec);

  std::cout << "cavity radius" << cavity->getRadius()[0] << std::endl;

}



void ReactionPotential::calc_rho_eff() {
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


}

void ReactionPotential::calc_gamma() {
  auto d_V_n = mrcpp::gradient(*derivative, V_n_func.real());
  gamma_func.alloc(NUMBER::Real);

  if(cavity->islinear() == true){
    println(0, "\nrunning linear");
    QMFunction temp_func;
    temp_func.alloc(NUMBER::Real);
    mrcpp::dot(this->apply_prec, temp_func.real(), d_V_n, d_cavity);
    qmfunction::multiply(gamma_func, temp_func, inv_eps_func, this->apply_prec);

  }else{
    println(0, "\nrunning exp.");
    mrcpp::dot(this->apply_prec, gamma_func.real(), d_V_n, d_cavity);

  }

  gamma_func.rescale(1.0/(4.0*MATHCONST::pi));

}

void ReactionPotential::setup(double prec) {
  setApplyPrec(prec);
  setup_eps();
  calc_rho_eff();

  d_cavity = mrcpp::gradient(*derivative, cavity_func.real());
  V_n_func.alloc(NUMBER::Real);
  mrcpp::apply(prec, V_n_func.real(), *poisson, rho_eff_func.real());
  double error = 1;
  int i = 1;
  while(error >= this->apply_prec) {
    gamma_func.free(NUMBER::Real);
    calc_gamma();
    QMFunction temp_func;
    qmfunction::add(temp_func, 1.0, rho_eff_func, 1.0, gamma_func, -1.0);
    QMFunction V_np1_func;
    V_np1_func.alloc(NUMBER::Real);
    mrcpp::apply(prec, V_np1_func.real(), *poisson, temp_func.real());

    QMFunction diff_func;
    qmfunction::add(diff_func, 1.0, V_n_func, -1.0, V_np1_func, -1.0);
    error = diff_func.norm();
    V_n_func = V_np1_func;

  //testing
    std::cout << "iterations: " << i << std::endl;
    println(0,"\nint. gamma \t" << gamma_func.integrate());
    i++;
  //testing
  }



  QMFunction V_0_func;


  V_0_func.alloc(NUMBER::Real);
  mrcpp::apply(prec, V_0_func.real(), *poisson, rho_tot.real());

  qmfunction::add(*this, 1.0, V_n_func, -1.0, V_0_func, -1.0);

  //testing

  QMFunction test_func;
  test_func.alloc(NUMBER::Real);
  qmfunction::add(test_func, 1.0, rho_eff_func, -1.0, rho_tot, -1.0);
  println(0, "\nint. rho_eff - rho\t" << test_func.integrate());

  for (double i = 0.1; i <= 10.00; i += 0.1) {
    double test1 = V_n_func.real().evalf({0.0, 0.0, i});
    double test2 = V_0_func.real().evalf({0.0, 0.0, i});
    double test3 = gamma_func.real().evalf({0.0, 0.0, i});
    double test4 = rho_eff_func.real().evalf({0.0, 0.0, i});
    double test5 = rho_tot.real().evalf({0.0, 0.0, i});
    double test6 = inv_eps_func.real().evalf({0.0, 0.0, i});
    double test7 = cavity->evalf({0.0, 0.0, i});

    std::cout << "func evalf at:\t" << i << "\t" << test1 << ' ' << test2 << ' ' << test3 << ' ' << test4 << ' ' << test5 << ' ' << test6 << " " << test7  << std::endl;

  }

  //testing

  V_n_func.free(NUMBER::Real);
  gamma_func.free(NUMBER::Real);

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

  QMFunction::free(NUMBER::Total);
}

} //namespace mrchem
