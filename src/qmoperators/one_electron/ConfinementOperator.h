#pragma once

#include "tensor/RankZeroOperator.h"
#include "ConfinementPotential.h"
#include "chemistry/Molecule.h"

/** @class ConfinementOperator
 * 
 * @brief Operator containing a ConfinementPotential
 *
 * A TensorOperator realization of @class ConfinementPotential.
 *
 */

namespace mrchem {

  class ConfinementOperator final : public RankZeroOperator {
  public:
    ConfinementOperator(double r_0, const int N, double s, std::vector<double> R, std::vector<mrcpp::Coord<3>> centers) {
      potential = std::make_shared<ConfinementPotential>(r_0, N, s, R, centers);

      // Invoke operator= to assign *this operator
      RankZeroOperator &Co = (*this);
      Co = potential;
      Co.name() = "Co";
    }
    ~ConfinementOperator() override = default;

    auto getRadius() { return this->potential->getRadius(); }
    auto getStiffness() { return this->potential->getStiffness(); }

  private:
    std::shared_ptr<ConfinementPotential> potential{nullptr};
  };

} // namespace mrchem
