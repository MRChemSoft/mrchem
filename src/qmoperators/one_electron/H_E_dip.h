#pragma once

#include "PositionOperator.h"
#include "chemistry/Nucleus.h"

/** Class H_E_dip 
 *
 * @brief Electric dipole operator
 *
 * The operator is based on the PositionOperator. The sign is coherent
 * with nuclear and electronic charges. It impleents explicit trace
 * functions for the nuclear contributions.
 *
 */

namespace mrchem {

/** @brief constructor
 *
 * @param[in] o the origin wrt which the dipole moment is computed.
 *
 */
class H_E_dip final : public PositionOperator {
public:
    H_E_dip(const mrcpp::Coord<3> &o)
            : PositionOperator(o) {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = -1.0 * (*this)[0];
        h[1] = -1.0 * (*this)[1];
        h[2] = -1.0 * (*this)[2];
    }

    /** @brief returns the nuclear contribution to the dipole moment
 *
 * @param[in] the set of nuclei
 *
 */
    ComplexVector trace(const Nuclei &nucs) {
        ComplexVector result = ComplexVector::Zero(3);
        for (auto &nuc_k : nucs) result += trace(nuc_k);
        return result;
    }

    /** @brief returns the contribution to the dipole moment
 *
 * @param[in] the nucleus
 *
 */
    ComplexVector trace(const Nucleus &nuc) {
        ComplexVector result = ComplexVector::Zero(3);
        double Z = nuc.getCharge();
        const mrcpp::Coord<3> &R = nuc.getCoord();
        result(0) = Z * this->r_x.function().real().evalf(R);
        result(1) = Z * this->r_y.function().real().evalf(R);
        result(2) = Z * this->r_z.function().real().evalf(R);
        return result;
    }

    using RankOneTensorOperator<3>::trace;
};

} //namespace mrchem
