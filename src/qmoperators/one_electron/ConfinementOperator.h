#pragma once

#include "tensor/RankZeroOperator.h"
#include "ConfinementPotential.h"

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
    ConfinementOperator(const double r_0, const int N) {
      potential = std::make_shared<ConfinementPotential>(r_0, N);

      // Invoke operator= to assign *this operator
      RankZeroOperator &Co = (*this);
      Co = potential;
      Co.name() = "Co";
    }
    ~ConfinementOperator() override = default;

    auto getRadius() { return this->potential->getRadius(); }
    auto getParam() { return this->potential->getParam(); }

  private:
    std::shared_ptr<ConfinementPotential> potential{nullptr};
  };

} // namespace mrchem
