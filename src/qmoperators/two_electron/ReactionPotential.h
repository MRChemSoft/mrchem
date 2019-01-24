#pragma once
#include "chemistry/Cavity.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "chemistry/chemistry_fwd.h"
#include "chemistry/Nucleus.h"

using namespace mrcpp;

namespace mrchem {


//start with a cavity initialized with a geometry and a standard gaussian rho with A*exp(-B*r^2) with A = (B/pi)^(3/2) include the orbitals later

class ReactionPotential final : public QMPotential{
public:

  ReactionPotential(mrcpp::PoissonOperator *P, mrcpp::DerivativeOperator<3> *D, Cavity *C, const Nuclei &nucs, OrbitalVector *Phi);
  ~ReactionPotential()= default;
  //void do_setup(double prec) { this->setup(prec); }

  double &get_tot_Energy();
  double &get_e_Energy();
  double &get_nuc_Energy();


  friend class ReactionOperator;

protected:
  void clear();


private:

  Cavity *cavity;
  Nuclei nuclei;
  OrbitalVector *orbitals;
  mrcpp::PoissonOperator *poisson;
  mrcpp::DerivativeOperator<3> *derivative;

  QMFunction rho_tot;
  QMFunction rho_el;
  QMFunction rho_nuc;
  QMFunction V_n_func;

  double e_Energy;
  double nuc_Energy;
  double tot_Energy;

  QMFunction setup_eps(bool is_eps);
  QMFunction calc_gamma(QMFunction inv_eps_func, mrcpp::FunctionTreeVector<3> d_cavity);
  QMFunction calc_rho_eff(QMFunction inv_eps_func);
  void setup(double prec);

};



} //namespace mrchem
