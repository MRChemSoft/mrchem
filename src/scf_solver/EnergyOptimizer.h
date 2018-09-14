#pragma once

#include "GroundStateSolver.h"

/** @class EnergyOptimizer
 *
 * @brief Ground state SCF optimization without kinetic energy operator
 *
 * This ground state SCF solver completely avoids the use of the kinetic energy
 * operator in the solution algorithm, in particular by computing a potential
 * update to the Fock matrix at the given iteration. This algorithm is considerably
 * less robust than the regular OrbitalOptimizer since there is no KAIN accelerator,
 * and is meant to be used only in the final iterations of the SCF procedure in order
 * to ensure quadratic precision of the energy relative to the orbitals (the effect
 * of this is currently questionable, after the improved implementation of derivative
 * operators in MRCPP).
 */

namespace mrchem {

class EnergyOptimizer final : public GroundStateSolver {
public:
    EnergyOptimizer(HelmholtzVector &h);

    void setup(FockOperator &fock,
               OrbitalVector &phi,
               ComplexMatrix &F,
               FockOperator &fock_np1,
               OrbitalVector &phi_np1);
    void clear();

    bool optimize();

protected:
    FockOperator *fOper_np1;     ///< Next iteration Fock operator (pointer to external object)
    OrbitalVector *orbitals_np1; ///< Next iteration orbitals (pointer to external object)

    ComplexMatrix calcFockMatrixUpdate(double prec, OrbitalVector &dPhi_n);
};

} // namespace mrchem
