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
    , V_n_func(false) {}

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
  qmfunction::add(rho_tot, 1.0, rho_e, 1.0, rho_N, -1.0);

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
  while(error >= this->apply_prec) {
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
  }
  QMFunction V_0_func;
  QMFunction V_eff_func;
  V_eff_func.alloc(NUMBER::Real);
  V_0_func.alloc(NUMBER::Real);
  mrcpp::apply(prec, V_0_func.real(), *poisson, rho_func.real());
  qmfunction::add(V_eff_func, 1.0, V_n_func, -1.0, V_0_func, -1.0);
  //V_eff = V_n_func - V_0_func
}
/*
void ReactionPotential::calc_d_Cavity(double prec){
  ABGVOperator &D = *this->*derivative;

  FunctionTree<3> dxC_tree(MRA);
  FunctionTree<3> dyC_tree(MRA);
  FunctionTree<3> dzC_tree(MRA);

  mrcpp::apply(dxC_tree, D, Cavity_tree, 0);
  mrcpp::apply(dyC_tree, D, Cavity_tree, 1);
  mrcpp::apply(dzC_tree, D, Cavity_tree, 2);

  d_cavity.pushback(std::make_tuple(1, &dxC_tree));
  d_cavity.pushback(std::make_tuple(1, &dyC_tree));
  d_cavity.pushback(std::make_tuple(1, &dzC_tree));

}

void ReactionPotential::calc_gamma(double prec){
  mrcpp::FunctionTreeVector<3> d_V;





}


  FunctionTree<3> *rho_eff_tree;
  FunctionTree<3> *gamma_tree;
  FunctionTree<3> *V_n_tree;
  FunctionTreeVector<3> *d_Cavity;*/

  void setup(double prec){}

/*  void calc_rho_eff(double prec);
  void calc_gamma(double prec);*/



} //namespace mrchem
