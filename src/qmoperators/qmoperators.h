#pragma once

#include "mrchem.h"

/* NOTES ON OPERATORS:
 *
 * The QMOperator is a fundamental quantum mechanical operator \hat{q} that can
 * be applied to an Orbital. However, its application is a private and can only
 * be done though a TensorOperator.
 *
 * The RankZeroTensorOperator is a general expansion of fundamental QMOperators:
 *      \hat{Q} = c_1 \hat{q}_{11}\hat{q}_{12}\cdots
 *              + c_2 \hat{q}_{21}\hat{q}_{22}\cdots
 *              + \cdots
 *
 * The RankOneTensorOperator is a vector of rank 0 operators:
 *      \hat{Q} = [Q_1, Q_2, \cdots, Q_N]
 *
 * The RankTwoTensorOperator is a matrix of rank 0 operators:
 *      \hat{Q} = [Q_11,   Q_12,   \cdots, Q_1N  ]
 *                [Q_21,   Q_22,   \cdots, Q_2N  ]
 *                [\cdots, \cdots, \cdots, \cdots]
 *                [Q_M1,   Q_M2,   \cdots, Q_MN  ]
 * 
 * A particular component can be applied in the following way:
 *      \psi = Q[i][j](\phi)
 *
 * This means that ALL operators needs to be implemented as a TensorOperator
 * containing at least one fundamental QMOperator (even those containing only
 * a single fundamental operator), but these can easily be reused by other
 * operators, e.g. the AngularMomentumOperator is constructed from
 * PositionOperator and MomentumOperator:
 * 
 *      l[0] = (r[1]p[2] - r[2]p[1])
 *      l[1] = (r[2]p[0] - r[0]p[2])
 *      l[2] = (r[0]p[1] - r[1]p[0])
 *
 * where position and momentum contain a single fundamental QMOperator for each
 * component (x,y,z).
 *
 */

namespace mrchem {

class QMOperator;
class RankZeroTensorOperator;
template<int I> class RankOneTensorOperator;
template<int I, int J> class RankTwoTensorOperator;

} //namespace mrchem
