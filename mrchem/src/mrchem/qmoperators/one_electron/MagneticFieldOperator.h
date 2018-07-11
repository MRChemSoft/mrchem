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

class MagneticFieldOperator final : public ExternalFieldOperator {
public:
 MagneticFieldOperator(const Eigen::Vector3d &f, mrcpp::DerivativeOperator<3> &D, const double *o = 0)
     : field(f), dipole(D, o) {
        RankZeroTensorOperator &d_x = this->dipole[0];
        RankZeroTensorOperator &d_y = this->dipole[1];
        RankZeroTensorOperator &d_z = this->dipole[2];

        RankZeroTensorOperator &HMF = (*this);
        HMF = f[0]*d_x + f[1]*d_y + f[2]*d_z;
    }
    
    ComplexDouble trace(const Nuclei &nucs) {
        return 0.0;
    }

    ComplexDouble trace(const Nucleus &nuc) {
        return 0.0;
    }

protected:
    Eigen::Vector3d field;
    H_B_dip dipole;
    
};

} //namespace mrchem
