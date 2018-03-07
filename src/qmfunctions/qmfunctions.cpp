#include "MRCPP/Printer"

#include "parallel.h"

#include "qmfunctions.h"
#include "Orbital.h"
#include "OrbitalVector.h"

using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace orbital {

/****************************************
 * Orbital related standalone functions *
 ****************************************/

/** Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
  * Complicated by the fact that both bra and ket can be interpreted as complex
  * conjugate versions of themselves. */
ComplexDouble dot(Orbital bra, Orbital ket) {
    if ((bra.spin() == SPIN::Alpha) and (ket.spin() == SPIN::Beta)) return 0.0;
    if ((bra.spin() == SPIN::Beta) and (ket.spin() == SPIN::Alpha)) return 0.0;

    double bra_conj(1.0), ket_conj(1.0);
    if (bra.conjugate()) bra_conj = -1.0;
    if (ket.conjugate()) ket_conj = -1.0;

    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra.hasReal() and ket.hasReal()) rr = mrcpp::dot(bra.real(), ket.real());
    if (bra.hasReal() and ket.hasImag()) ri = mrcpp::dot(bra.real(), ket.imag());
    if (bra.hasImag() and ket.hasReal()) ir = mrcpp::dot(bra.imag(), ket.real());
    if (bra.hasImag() and ket.hasImag()) ii = mrcpp::dot(bra.imag(), ket.imag());

    double real_part = rr + bra_conj*ket_conj*ii;
    double imag_part = ket_conj*ri - bra_conj*ir;
    return ComplexDouble(real_part, imag_part);
}

/** out = a*inp_a + b*inp_b
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves. */
Orbital add(ComplexDouble a, Orbital inp_a, ComplexDouble b, Orbital inp_b, double prec) {
    ComplexVector coefs(2);
    coefs(0) = a;
    coefs(1) = b;

    OrbitalVector orbs;
    orbs.push_back(inp_a);
    orbs.push_back(inp_b);

    return add(coefs, orbs, prec);
}

/** out_i = a*(inp_a)_i + b*(inp_b)_i
 *
 *  Component-wise addition of orbitals.
 *
 */
OrbitalVector add(ComplexDouble a, OrbitalVector inp_a,
                  ComplexDouble b, OrbitalVector inp_b,
                  double prec) {
    if (inp_a.size() != inp_b.size()) MSG_ERROR("Size mismatch");

    OrbitalVector out;
    for (int i = 0; i < inp_a.size(); i++) {
        Orbital out_i = add(a, inp_a[i], b, inp_b[i]);
        out.push_back(out_i);
    }
    return out;
}

/** out = c_0*inp_0 + c_1*inp_1 + ...
  *
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves.
  *
  */
Orbital add(const ComplexVector &c, OrbitalVector inp, double prec) {
    if (c.size() != inp.size()) MSG_ERROR("Size mismatch");

    Orbital out;
    // set output spin from first contributing input
    for (int i = 0; i < inp.size(); i++) {
        if (abs(c[i]) < mrcpp::MachineZero) continue;
        out = inp[i].paramCopy();
        break;
    }
    // all contributing input spins must be equal
    for (int i = 0; i < inp.size(); i++) {
        if (abs(c[i]) < mrcpp::MachineZero) continue;
        if (out.spin() != inp[i].spin()) MSG_FATAL("Mixing spins");
    }

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    double thrs = mrcpp::MachineZero;
    for (int i = 0; i < inp.size(); i++) {
        bool cHasReal = (abs(c[i].real()) > thrs);
        bool cHasImag = (abs(c[i].imag()) > thrs);

        double conj(1.0);
        if (inp[i].conjugate()) conj = -1.0;

        if (cHasReal and inp[i].hasReal()) rvec.push_back(      c[i].real(), &inp[i].real());
        if (cHasImag and inp[i].hasImag()) rvec.push_back(-conj*c[i].imag(), &inp[i].imag());

        if (cHasImag and inp[i].hasReal()) ivec.push_back(      c[i].imag(), &inp[i].real());
        if (cHasReal and inp[i].hasImag()) ivec.push_back( conj*c[i].real(), &inp[i].imag());
    }

    if (rvec.size() > 0) {
        out.alloc(NUMBER::Real);
        if (prec < 0.0) {
            mrcpp::copy_grid(out.real(), rvec);
            mrcpp::add(prec, out.real(), rvec, 0);
        } else {
            mrcpp::add(prec, out.real(), rvec);
        }
    }
    if (ivec.size() > 0) {
        out.alloc(NUMBER::Imag);
        if (prec < 0.0) {
            mrcpp::copy_grid(out.imag(), ivec);
            mrcpp::add(prec, out.imag(), ivec, 0);
        } else {
            mrcpp::add(prec, out.imag(), ivec);
        }
    }
    return out;
}

/* Orbital transformation out_vec = U*inp_vec
 *
 * The transformation matrix is not necessarily square.
 *
 */
OrbitalVector add(const ComplexMatrix &U, OrbitalVector inp, double prec) {
    if (mpi::orb_size > 1) NOT_IMPLEMENTED_ABORT;

    OrbitalVector out;
    for (int i = 0; i < U.rows(); i++) {
        const ComplexVector &c = U.row(i);
	Orbital out_i = add(c, inp, prec);
        out.push_back(out_i);
    }
    return out;
}


/** out = inp_a * inp_b
  *
  * Complicated by the fact that both inputs can be interpreted
  * as complex conjugate versions of themselves.
  *
  */
Orbital multiply(Orbital inp_a, Orbital inp_b, double prec) {
    int occ = compare_occ(inp_a, inp_b);
    int spin = compare_spin(inp_a, inp_b);
    Orbital out(spin, occ);

    double a_conj(1.0), b_conj(1.0);
    if (inp_a.conjugate()) a_conj = -1.0;
    if (inp_b.conjugate()) b_conj = -1.0;

    { // Real part
        FunctionTreeVector<3> vec;
        if (inp_a.hasReal() and inp_b.hasReal()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = 1.0;
            if (prec < 0.0) {
                // Union grid
                mrcpp::copy_grid(*tree, inp_a.real());
                mrcpp::copy_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real());
            }
            vec.push_back(tree);
        }
        if (inp_a.hasImag() and inp_b.hasImag()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = -1.0*a_conj*b_conj;
            if (prec < 0.0) {
                mrcpp::copy_grid(*tree, inp_a.imag());
                mrcpp::copy_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag(), 0);
            } else {
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag());
            }
            vec.push_back(tree);
        }
        if (vec.size() == 1) {
            out.setReal(vec[0]);
            vec.clear(false);
        }
        if (vec.size() == 2) {
            out.alloc(NUMBER::Real);
            mrcpp::copy_grid(out.real(), vec);
            mrcpp::add(prec, out.real(), vec, 0);
            vec.clear(true);
        }
    }

    { // Imaginary part
        FunctionTreeVector<3> vec;
        if (inp_a.hasReal() and inp_b.hasImag()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = b_conj;
            if (prec < 0.0) {
                // Union grid
                mrcpp::copy_grid(*tree, inp_a.real());
                mrcpp::copy_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag());
            }
            vec.push_back(tree);
        }
        if (inp_a.hasImag() and inp_b.hasReal()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = a_conj;
            if (prec < 0.0) {
                // Union grid
                mrcpp::copy_grid(*tree, inp_a.imag());
                mrcpp::copy_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real());
            }
            vec.push_back(tree);
        }
        if (vec.size() == 1) {
            out.setImag(vec[0]);
            vec.clear(false);
        }
        if (vec.size() == 2) {
            out.alloc(NUMBER::Imag);
            mrcpp::copy_grid(out.imag(), vec);
            mrcpp::add(prec, out.imag(), vec, 0);
            vec.clear(true);
        }
    }

    return out;
}

/** Compare spin and occupancy of two orbitals
 *
 *  Returns true if orbital parameters are the same.
 *
 */
bool compare(const Orbital &orb_a, const Orbital &orb_b) {
    bool comp = true;
    if (compare_occ(orb_a, orb_b) < 0) {
        MSG_WARN("Different occupancy");
        comp = false;
    }
    if (compare_spin(orb_a, orb_b) < 0) {
        MSG_WARN("Different spin");
        comp = false;
    }
    return comp;
}

/** Compare occupancy of two orbitals
 *
 *  Returns the common occupancy if they match, -1 if they differ.
 *
 */
int compare_occ(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.occ() == orb_b.occ()) comp = orb_a.occ();
    return comp;
}

/** Compare spin of two orbitals
 *
 *  Returns the common spin if they match, -1 if they differ.
 *
 */
int compare_spin(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.spin() == orb_b.spin()) comp = orb_a.spin();
    return comp;
}

} //namespace orbital


namespace density {

/****************************************
 * Density related standalone functions *
 ****************************************/

} //namespace density

} //namespace mrchem
