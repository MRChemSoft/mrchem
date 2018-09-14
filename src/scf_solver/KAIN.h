#pragma once

#include "Accelerator.h"

/** @class KAIN
 *
 * This class implements the Krylov subspace Accelerated Inexact Newton (KAIN)
 * method as described by R.J. Harrison (J. Comput. Chem. 25, 328, 2004).
 *
 * The SCF problem is formulated as finding the root of the function
 *
 * \f$ f(x^n) = -2H[x^n] - x^n \f$
 *
 * where \f$ x^n \f$ is the vector of orbitals, possibly appended by the
 * Fock matrix \f$ x^n = (\phi^n_0, \phi^n_1, \dots \phi^n_N, F^n) \f$
 */

namespace mrchem {

class KAIN final : public Accelerator {
public:
    KAIN(int max, int min = 0, bool sep = false)
            : Accelerator(max, min, sep) {}

protected:
    void setupLinearSystem();
    void expandSolution(double prec, OrbitalVector &Phi, OrbitalVector &dPhi, ComplexMatrix *F, ComplexMatrix *dF);
};

} // namespace mrchem
