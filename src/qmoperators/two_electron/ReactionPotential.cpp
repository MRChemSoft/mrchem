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

ReactionPotential::ReactionPotential( mrcpp::PoissonOperator *P, mrcpp::DerivativeOperator<3> *D, Cavity *C, const Nuclei &nucs, OrbitalVector *Phi)
    : QMPotential(1, false)
    , cavity(C)
    , nuclei(nucs)
    , orbitals(Phi)
    , poisson(P)
    , derivative(D)
    , rho_func(false)
    , cavity_func(false)
    , inv_eps_func(false)
    , rho_eff_func(false)
    , gamma_func(false)
    , V_n_func(false)
    , V_eff_func(false) {}

//  ~ReactionPotential();

void ReactionPotential::setup_eps() {

  cavity->eval_epsilon(false, cavity->is_linear);

  qmfunction::project(cavity_func, *cavity, NUMBER::Real, this->apply_prec);
  cavity_func.rescale(cavity->dcoeff);

  cavity->eval_epsilon(true, cavity->is_linear);

  qmfunction::project(inv_eps_func, *cavity, NUMBER::Real, this->apply_prec);

}



void ReactionPotential::calc_rho_eff() {
  Density rho_N = chemistry::compute_nuclear_density(this->apply_prec, this->nuclei, 1000);
  Density rho_e(false);
  Density rho_tot(false);
  density::compute(this->apply_prec, rho_e, *orbitals, DENSITY::Total);
  qmfunction::add(rho_tot, -1.0, rho_e, 1.0, rho_N, -1.0);

  this->rho_func = rho_tot;

  qmfunction::multiply(rho_eff_func, rho_tot, inv_eps_func, this->apply_prec);


}

void ReactionPotential::calc_gamma() {
  auto d_V_n = mrcpp::gradient(*derivative, V_n_func.real());
  gamma_func.alloc(NUMBER::Real);
  mrcpp::dot(this->apply_prec, gamma_func.real(), d_V_n, d_cavity);
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
  while(error >= 10*this->apply_prec) {
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
    if(error <= this->apply_prec){
      println(0, "\nR_char:\t" << gamma_func.integrate());
    }
    gamma_func.free(NUMBER::Real);
    i++;
  }

  QMFunction V_0_func;


  V_0_func.alloc(NUMBER::Real);
  mrcpp::apply(prec, V_0_func.real(), *poisson, rho_func.real());

  QMFunction temp_prod_func;
  qmfunction::add(temp_prod_func, 1.0, V_n_func, -1.0, V_0_func, -1.0);

  qmfunction::multiply(V_eff_func, rho_func, temp_prod_func, this->apply_prec);

  println(0, -0.5/2.0);

  println(0, "\ncharge:\t" << rho_func.integrate());
  println(0, "E_r:\t" << V_eff_func.integrate());
}

double ReactionPotential::getEnergy(){
 return 0.5*this->V_eff_func.integrate().real();
}

void ReactionPotential::clear() {
  clearApplyPrec();
  //QMFunction::free(NUMBER::Total);



}

} //namespace mrchem
