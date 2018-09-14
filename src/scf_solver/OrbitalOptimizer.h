#pragma once

#include "GroundStateSolver.h"

/** @class OrbitalOptimizer
 *
 * @brief Ground state SCF optimization with kinetic energy operator
 *
 * This ground state SCF solver computes the Fock matrix explicitly by application of
 * the kinetic energy operator. This simplifies the algorithm significanly when a KAIN
 * iterative accelerator is present, or when diagonalization/localization occurs.
 * The derivative operator in the kinetic energy MIGHT affect the quadratic precision
 * in the energy, but this point seems to no longer be critical.
 */

namespace mrchem {

class Accelerator;

class OrbitalOptimizer final : public GroundStateSolver {
public:
    OrbitalOptimizer(HelmholtzVector &h, Accelerator *k = 0);

    void setup(FockOperator &fock, OrbitalVector &phi, ComplexMatrix &F);
    void clear();

    bool optimize();

protected:
    Accelerator *kain; ///< KAIN accelerator(pointer to external object)
};

} // namespace mrchem
