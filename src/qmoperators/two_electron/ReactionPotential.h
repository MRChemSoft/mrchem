#include "chemistry/Cavity.h"
#include "qmoperators/one_electron/QMPotential.h"


using namespace mrcpp;

namespace mrchem{


//start with a cavity initialized with a geometry and a standard gaussian rho with A*exp(-B*r^2) with A = (B/pi)^(3/2) include the orbitals later

class ReactionPotential final : public QMPotential{
public:
  ReactionPotential(Cavity *cav,  mrcpp::PoissonOperator *P, mrcpp::ABGVOperator<3> *D);
  ~ReactionPotential();


protected:
  Cavity *Cavity;  

  mrcpp::PoissonOperator *poisson;
  mrcpp::ABGVOperator<3> *derivative;
  
  mrcpp::FunctionTree<3> *Cavity_tree;
  mrcpp::FunctionTree<3> *inv_eps_tree;
  mrcpp::FunctionTree<3> *rho_eff_tree;
  mrcpp::FunctionTree<3> *gamma_tree;
  mrcpp::FunctionTree<3> *V_n_tree;
  mrcpp::FunctionTreeVector<3> *d_Cavity;
  
  void setup(double prec);
  
  void setup_eps(double prec);
  void calc_rho_eff(double prec);
  void calc_gamma(double prec);
};



} //namespace mrchem
