#pragma once

#include "RankZeroTensorOperator.h"
#include "H_E_dip.h"

/** @class ElectricFieldOperator
 *
 * @brief External electric field operator
 *
 * An external electric field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the electric field vector.
 *
 */

namespace mrchem {

class ElectricFieldOperator final : public RankZeroTensorOperator {
public:
 ElectricFieldOperator(const Eigen::Vector3d &f)
     : field(f) {
        RankZeroTensorOperator &d_x = this->dipole[0];
        RankZeroTensorOperator &d_y = this->dipole[1];
        RankZeroTensorOperator &d_z = this->dipole[2];

        RankZeroTensorOperator &HEF = (*this);
        HEF = + f[0]*d_x + f[1]*d_y + f[2]*d_z;
    }

    ComplexDouble trace(const Nuclei &nucs) {
        ComplexDouble result = 0.0;
        for (int k = 0; k < nucs.size(); k++) {
            result += trace(nucs[k]);
        }
        return result;
    }

    ComplexDouble trace(const Nucleus &nuc) {
        return dipole.trace(nuc).dot(field);
    }

    using RankZeroTensorOperator::trace;

 protected:
    Eigen::Vector3d field;
    H_E_dip dipole;
    

};

} //namespace mrchem
