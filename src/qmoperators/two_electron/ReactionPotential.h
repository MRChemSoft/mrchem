#include "chemistry/Cavity.h"
#include "FunctionTree.h"
#include <vector>
#include <array>
using namespace mrcpp;

namespace mrchem{

class ReactionPotential{
public:
  ReactionPotential(const FunctionTree<3> &rho_tree, FunctionTree<3> cavity_tree, FunctionTree<3> epsilon_tree); //initialize in default
  void rho_eff();
  void gamma();
  
protected:
  FunctionTree<3> rho_tree, cavity_tree, epsilon_tree, rho_eff_tree, gamma_tree;
  bool is_lin = false;
  
  


};






} //namespace mrchem
