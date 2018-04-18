#include "utils/NonlinearMaximizer.h"

#include "qmfunctions.h"

/** subclass which defines the particular Gradient and Hessian
 * and other specific functions for a maximization of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 *
 */

namespace mrchem {

class RRMaximizer final : public NonlinearMaximizer {
public:
    RRMaximizer(double prec, OrbitalVector &Phi);
    const DoubleMatrix &getTotalU() const { return this->total_U; }

protected:
    int N;                  // number of orbitals
    DoubleMatrix r_i_orig;  // <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
    DoubleMatrix r_i ;      // rotated  r_i_orig
    DoubleMatrix total_U;   // the rotation matrix of the orbitals

    double functional() const;
    double make_gradient();
    double make_hessian();
    void do_step(const DoubleVector &step);
};

} //namespace mrchem
