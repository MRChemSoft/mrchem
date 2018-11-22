#include <cmath>
#include "MRCPP/MWOperators"
#include "MRCPP/Gaussians"

#include "Cavity.h"
#include "ReactionPotential.h"

using mrcpp::FunctionTree;
using mrcpp::PoissonOperator;
using mrcpp::ABGVOperator;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

ReactionPotential::ReactionPotential(Cavity *cav,  PoissonOperator *P, ABGVOperator<3> *D)
  	: Cavity(cav),
	  poisson(P),
	  derivative(D) {}

//  ~ReactionPotential();

void ReactionPotential::setup_eps(double prec){
  Cavity.eval_epsilon(false, Cavity.is_linear);
  mrcpp::project(prec, Cavity_tree, Cavity.evalf);
  Cavity_tree.rescale(Cavity.dcoeff);

  Cavity.eval_epsilon(true, Cavity.is_linear);
  mrcpp::project(prec, inv_eps_tree, Cavity.evalf);

}

void ReactionPotential::calc_rho_eff(double prec) {
  
  // make rho as a gaussian
  double alpha, beta;
  const double pi = std::atan(1.0)*4;
  beta = 100;
  alpha = beta/pi * std::sqrt(beta/pi);
  double pos[3] = {0, 0, 0};
  int pow[3] = {0, 0, 0};
  mrcpp::GaussFunc<3> rho(beta, alpha, pos, pow);
  FunctionTree<3> rho_tree(MRA);
  mrcpp::build_grid(rho_tree, rho);
  mrcpp::project(prec, rho_tree, rho);
  
  // might lose accuracy
  // calculate rho_eff
  mrcpp::multiply(prec, rho_eff_tree, 1, rho_tree, inv_eps_tree);
  }

}

void ReactionPotential::calc_d_Cavity(double prec){
  ABGVOperator &D = *this->derivative;

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
  FunctionTreeVector<3> *d_Cavity;

  void setup(double prec);

  void calc_rho_eff(double prec);
  void calc_gamma(double prec);



} //namespace mrchem






