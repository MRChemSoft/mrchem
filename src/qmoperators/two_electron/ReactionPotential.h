#include "chemistry/Cavity.h"
#include "FunctionTree.h"
#include <vector>
#include <array>
using namespace mrcpp;

namespace mrchem{

class ReactionPotential{
public:
  ReactionPotential(const FunctionTree<3> &rho_tree, std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope); //initialize in default
  ReactionPotential(const FunctionTree<3> &rho_tree, std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope, double e_i, double e_o, bool is_lin);
  void rho_eff();
  void gamma();

protected:
  bool is_lin = false;
  Cavity C;
  Cavity eps_inv;
  FunctionTree<3> rho_tree;
  FunctionTree<3> rho_eff_tree;
  FunctionTree<3> gamma_tree;

  double prec = 1.0e-6;
  int n = -4;
  int l[3] = {-1, -1, -1};
  int nb[3] = {2, 2, 2};
  int N = 25;
  int k = 7;

  BoundingBox<3> world(n, l, nb);
  ScalingBasis basis =  InterpolatingBasis(k);
  MultiResolutionAnalysis<3> MRA(world, basis, N);


};






} //namespace mrchem
