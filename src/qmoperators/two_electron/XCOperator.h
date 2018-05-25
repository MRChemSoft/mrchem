#pragma once

#include "RankZeroTensorOperator.h"
#include "XCPotential.h"

/** @class XCOperator
 *
 * @brief DFT Exchange-Correlation operator containing a single XCPotential
 *
 * This class is a simple TensorOperator realization of @class XCPotential.
 *
 */

namespace mrchem {
    
class XCOperator final : public RankZeroTensorOperator {
public:
    XCOperator(mrdft::XCFunctional *F, OrbitalVector *Phi = nullptr)
            : potential(F, Phi) {
        RankZeroTensorOperator &XC = (*this);
        XC = this->potential;
    }

    double getEnergy() { return this->potential.getEnergy(); }
    Density &getDensity(int spin) { return this->potential.getDensity(spin); }

protected:
    XCPotential potential;
};

} //namespace mrchem
