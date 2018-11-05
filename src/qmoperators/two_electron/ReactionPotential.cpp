#include "chemistry/Cavity.h"
#include "ReactionPotentila.cpp"
using namespace mrcpp;

namespace mrchem {

ReactionPotential::ReactionPotential(const FunctionTree<3> &rho_tree, std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope){

  Cavity C(coord, R, slope);
  Cavity eps_inv(coord, R, slope);
  eps_inv.eval_epsilon(true, false);
  
  this->C = C;
  this->rho_tree = rho_tree;
  this->eps_inv = eps_inv;
}


ReactionPotential::ReactionPotential(const FunctionTree<3> &rho_tree, std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope, double e_i, double e_o, bool is_lin){
  Cavity C(coord, R, slope, e_i, e_o);
  Cavity eps_inv(coord, R, slope);
  eps_inv.eval_epsilon(true, is_lin);

  this->C = C;
  this->rho_tree = rho_tree;
  this->eps_inv = eps_inv;
  this->is_lin = is_lin;


}


void ReactionPotential::rho_eff(){
  FunctionTree<3> rho_eff_tree (MRA);
  FunctionTree<3> eps_inv_tree (MRA);
  build_grid(eps_inv_tree, eps_inv.evalf);
  project(prec, eps_inv_tree, eps_inv.evalf);
  
  multiply(prec, rho_eff_tree, 1.0, rho_tree, eps_inv_tree);
  this->rho_eff_tree = rho_eff_tree;
}


void ReactionPotential::gamma(){
  if(is_lin == false){
  

  }
}



} //namespace mrcem
