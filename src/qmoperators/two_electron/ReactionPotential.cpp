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

  cavity->eval_epsilon(false, cavity->islinear());

  qmfunction::project(cavity_func, *cavity, NUMBER::Real, this->apply_prec);
  cavity_func.rescale(cavity->dcoeff);

  cavity->eval_epsilon(true, cavity->islinear());

  qmfunction::project(inv_eps_func, *cavity, NUMBER::Real, this->apply_prec);

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

    QMFunction temp_func;
    temp_func.alloc(NUMBER::Real);
    mrcpp::dot(this->apply_prec, temp_func.real(), d_V_n, d_cavity);
    qmfunction::multiply(gamma_func, temp_func, inv_eps_func, this->apply_prec);

  }else{

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

    if(error >= 10*this->apply_prec) break;
    gamma_func.free(NUMBER::Real);
  }

  QMFunction V_0_func;


  V_0_func.alloc(NUMBER::Real);
  mrcpp::apply(prec, V_0_func.real(), *poisson, rho_tot.real());

  qmfunction::add(*this, 1.0, V_n_func, -1.0, V_0_func, -1.0);
  /*
  //testing
  println(0, "\nis it eps\t" << cavity->iseps());

  println(0, "\ncav coordinates\n")
  println(0, "atom number:\t" << 0);
  double tempx = cavity->getcoords()[0][0];
  double tempy = cavity->getcoords()[0][1];
  double tempz = cavity->getcoords()[0][2];
  println(0, tempx << ' ' << tempy << ' ' << tempz << "\n");

  println(0, "atom number:\t" << 1);
  double tempx = cavity->getcoords()[1][0];
  double tempy = cavity->getcoords()[1][1];
  double tempz = cavity->getcoords()[1][2];
  println(0, tempx << ' ' << tempy << ' ' << tempz << "\n");

  println(0, "\ncav at (0.0, 0.0, 0.0)\t" << cavity->evalf({0.0, 0.0, 0.0})); //gabriel: debugging
  println(0, "\ncav at (0.0, 0.0, 10.0)\t" << cavity->evalf({0.0, 0.0, 10.0})); //gabriel: debugging
  println(0, "\nE_r/2 should be:\t" << ((1.0 - 2.0)*pow(rho_el.integrate().real(), 2))/(4.0*2.0));
  println(0, "\nE_r/2 is:\t" << get_tot_Energy());*/

  V_n_func.free(NUMBER::Real);
  gamma_func.free(NUMBER::Real);

}

double &ReactionPotential::get_tot_Energy(){
  QMFunction temp_prod_func;
  qmfunction::multiply(temp_prod_func, rho_tot, *this, this->apply_prec);
  tot_Energy = temp_prod_func.integrate().real()/2;
  return tot_Energy;
}


double &ReactionPotential::get_e_Energy(){
  QMFunction temp_prod_func;
  qmfunction::multiply(temp_prod_func, rho_el, *this, this->apply_prec);
  e_Energy = temp_prod_func.integrate().real()/2;
  return e_Energy;
}


double &ReactionPotential::get_nuc_Energy(){
  QMFunction temp_prod_func;
  qmfunction::multiply(temp_prod_func, rho_nuc, *this, this->apply_prec);
  nuc_Energy = temp_prod_func.integrate().real()/2;
  return nuc_Energy;
}

void ReactionPotential::clear() {
  clearApplyPrec();

  QMFunction::free(NUMBER::Total);  



}

} //namespace mrchem





}

} //namespace mrchem
