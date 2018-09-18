#pragma once

#include <deque>
#include <vector>

#include "mrchem.h"

#include "qmfunctions/qmfunction_fwd.h"

/** @class Acccelerator
 *
 *  Base class for iterative subspace accelerators for use in SCF
 *  optimizations. Solves a linear system of equations \f$ Ac = b \f$
 *  to obtain the coefficient vector \f$ c \f$ that gives the linear
 *  combination of history orbitals that gives the next iteration.
 *  The size of the problem (matrix A) is given by the length of the
 *  history (not the number of orbitals).
 *
 *  The form of the linear problem (how to compute A and b) is given
 *  in the subclass. It is possible to treat each orbital separately,
 *  but usually you will benefit from including all orbitals in the
 *  same subspace and solve for them simultaneously (default). It is
 *  also possible to include the Fock matrix and corresponding update
 *  to the subspace. In this case one entry is added (corresponding to
 *  an extra orbital) using the Frobenius inner product of matrices.
 */

namespace mrchem {

class Accelerator {
public:
    Accelerator(int max, int min, bool sep);
    virtual ~Accelerator() { this->clear(); }
    void clear();

    void setMaxHistory(int max) { this->maxHistory = max; }
    void setMinHistory(int min) { this->minHistory = min; }

    // clang-format off
    void accelerate(double prec,
                    OrbitalVector &Phi,
                    OrbitalVector &dPhi,
                    ComplexMatrix *F = nullptr,
                    ComplexMatrix *dF = nullptr);
    // clang-format on

    void copyOrbitals(OrbitalVector &Phi, int nHistory = 0);
    void copyOrbitalUpdates(OrbitalVector &dPhi, int nHistory = 0);

    void replaceOrbitals(OrbitalVector &Phi, int nHistory = 0);
    void replaceOrbitalUpdates(OrbitalVector &dPhi, int nHistory = 0);

    void rotate(const ComplexMatrix &U, bool all = true);

    int printTreeSizes() const;

protected:
    int minHistory;   ///< Accelerator is activated when history reaches this size
    int maxHistory;   ///< Oldest iteration is discarded when history exceeds this size
    bool sepOrbitals; ///< Use separate subspace for each orbital

    std::vector<DoubleMatrix> A; ///< Vector of A matrices
    std::vector<DoubleVector> b; ///< Vector of b vectors
    std::vector<DoubleVector> c; ///< Vector of c vectors

    std::deque<OrbitalVector> orbitals;  ///< Orbital history
    std::deque<OrbitalVector> dOrbitals; ///< Orbital update history
    std::deque<ComplexMatrix> fock;      ///< Fock history
    std::deque<ComplexMatrix> dFock;     ///< Fock update history

    bool verifyOverlap(OrbitalVector &phi);

    // clang-format off
    void push_back(OrbitalVector &phi,
                   OrbitalVector &dPhi,
                   ComplexMatrix *F = nullptr,
                   ComplexMatrix *dF = nullptr);
    // clang-format on

    void solveLinearSystem();
    void clearLinearSystem();
    // clang-format off
    void sortLinearSystem(std::vector<DoubleMatrix> &A_mat,
                          std::vector<DoubleVector> &b_vec);
    // clang-format on

    virtual void setupLinearSystem() = 0;
    virtual void expandSolution(double prec,
                                OrbitalVector &phi,
                                OrbitalVector &dPhi,
                                ComplexMatrix *F,
                                ComplexMatrix *dF) = 0;
};

} // namespace mrchem
