#pragma once
#include "chemistry/Cavity.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "chemistry/chemistry_fwd.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/Density.h"

using namespace mrcpp;

namespace mrchem {


//start with a cavity initialized with a geometry and a standard gaussian rho with A*exp(-B*r^2) with A = (B/pi)^(3/2) include the orbitals later

class ReactionPotential final : public QMPotential{
public:

  ReactionPotential(mrcpp::PoissonOperator *P, mrcpp::DerivativeOperator<3> *D, Cavity *C, const Nuclei &nucs, OrbitalVector *Phi);
  ~ReactionPotential()= default;

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

  Density rho_tot;
  Density rho_el;
  Density rho_nuc;
  QMFunction V_n_func;

  double e_Energy;
  double nuc_Energy;
  double tot_Energy;

  void calc_eps(bool is_inv, QMFunction &cavity_func);
  void calc_rho_eff(QMFunction const &inv_eps_func, QMFunction &rho_eff_func);
  void calc_gamma(QMFunction const &inv_eps_func, QMFunction &gamma_func, mrcpp::FunctionTreeVector<3> &d_cavity);
  void grad_G(QMFunction &gamma_func, QMFunction &cavity_func, QMFunction &rho_tot, QMFunction &grad_G_func);
  void setup(double prec);

};



} //namespace mrchem
