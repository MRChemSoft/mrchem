#pragma once

#include "RankZeroTensorOperator.h"
#include "H_B_dip.h"

/** @class MagneticFieldOperator
 *
 * @brief External magnetic field operator
 *
 * An external magnetic field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the magnetic field vector.
 *
 */

namespace mrchem {

class MagneticFieldOperator final : public RankZeroTensorOperator {
public:
 MagneticFieldOperator(const Eigen::Vector3d &f, mrcpp::DerivativeOperator<3> &D, const double *o = 0)
     : field(f), dipole(D, o) {
        RankZeroTensorOperator &d_x = this->dipole[0];
        RankZeroTensorOperator &d_y = this->dipole[1];
        RankZeroTensorOperator &d_z = this->dipole[2];

        RankZeroTensorOperator &HMF = (*this);
        HMF = 0.5*(f[0]*d_x + f[1]*d_y + f[2]*d_z);
    }

protected:
    Eigen::Vector3d field;
    H_B_dip dipole;
    
};

} //namespace mrchem
