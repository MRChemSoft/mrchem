#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#pragma GCC system_header
#include <Eigen/Dense>

#include <deque>
#include <vector>

#include "OrbitalAdder.h"

/** Base class for iterative subspace accelerators for use in SCF
  * optimizations. Solves a linear system of equations \f$ Ac = b \f$
  * to obtain the coefficient vector \f$ c \f$ that gives the linear
  * combination of history orbitals that gives the next iteration.
  * The size of the problem (matrix A) is given by the length of the
  * history (not the number of orbitals).
  *
  * The form of the linear problem (how to compute A and b) is given
  * in the subclass. It is possible to treat each orbital separately,
  * but usually you will benefit from including all orbitals in the
  * same subspace and solve for them simultaneously (default). It is
  * also possible to include the Fock matrix and corresponding update
  * to the subspace. In this case one entry is added (corresponding to
  * an extra orbital) using the Frobenius inner product of matrices. */
class Accelerator {
public:
    Accelerator(int max, int min, bool sep);
    virtual ~Accelerator() { this->clear(); }
    void clear();

    void setMaxHistory(int max) { this->maxHistory = max; }
    void setMinHistory(int min) { this->minHistory = min; }

    void accelerate(double prec,
                    OrbitalVector &phi,
                    OrbitalVector &dPhi,
                    Eigen::MatrixXd *F = 0,
                    Eigen::MatrixXd *dF = 0);

    void copyOrbitals(OrbitalVector &phi, int nHistory = 0);
    void copyOrbitalUpdates(OrbitalVector &dPhi, int nHistory = 0);

    void replaceOrbitals(OrbitalVector &phi, int nHistory = 0);
    void replaceOrbitalUpdates(OrbitalVector &dPhi, int nHistory = 0);

    void rotate(const Eigen::MatrixXd &U, bool rotAll = true);

protected:
    int minHistory;   ///< Accelerator is activated when history reaches this size
    int maxHistory;   ///< Oldest iteration is discarded when history exceeds this size
    bool sepOrbitals; ///< Use separate subspace for each orbital
    OrbitalAdder add;

    std::vector<Eigen::MatrixXd *> A;   ///< Vector of A matrices
    std::vector<Eigen::VectorXd *> b;   ///< Vector of b vectors
    std::vector<Eigen::VectorXd *> c;   ///< Vector of c vectors

    std::deque<OrbitalVector *> orbitals;	 ///< Orbital history
    std::deque<OrbitalVector *> dOrbitals; ///< Orbital update history
    std::deque<Eigen::MatrixXd> fock;	 ///< Fock history
    std::deque<Eigen::MatrixXd> dFock;  ///< Fock update history

    bool verifyOverlap(OrbitalVector &phi);

    void push_back(OrbitalVector &phi,
                   OrbitalVector &dPhi,
                   Eigen::MatrixXd *F = 0,
                   Eigen::MatrixXd *dF = 0);

    void solveLinearSystem();
    void clearLinearSystem();
    void sortLinearSystem(std::vector<Eigen::MatrixXd *> &A_mat,
                          std::vector<Eigen::VectorXd *> &b_vec);

    virtual void setupLinearSystem() = 0;
    virtual void expandSolution(OrbitalVector &phi,
                                OrbitalVector &dPhi,
                                Eigen::MatrixXd *F,
                                Eigen::MatrixXd *dF) = 0;

};

