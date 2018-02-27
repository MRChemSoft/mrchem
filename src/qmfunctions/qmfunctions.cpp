#include "MRCPP/Printer"

#include "qmfunctions.h"
#include "Orbital.h"

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
    int occ = compare_occ(inp_a, inp_b);
    int spin = compare_spin(inp_a, inp_b);
    Orbital out(occ, spin);

    double thrs = mrcpp::MachineZero;
    bool aHasReal = (abs(a.real()) > thrs);
    bool aHasImag = (abs(a.imag()) > thrs);
    bool bHasReal = (abs(b.real()) > thrs);
    bool bHasImag = (abs(b.imag()) > thrs);

    double a_conj(1.0), b_conj(1.0);
    if (inp_a.conjugate()) a_conj = -1.0;
    if (inp_b.conjugate()) b_conj = -1.0;

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    if (inp_a.hasReal() and aHasReal) rvec.push_back(a.real(), &inp_a.real());
    if (inp_b.hasReal() and bHasReal) rvec.push_back(b.real(), &inp_b.real());
    if (inp_a.hasImag() and aHasImag) rvec.push_back(-a_conj*a.imag(), &inp_a.imag());
    if (inp_b.hasImag() and bHasImag) rvec.push_back(-b_conj*b.imag(), &inp_b.imag());

    if (inp_a.hasReal() and aHasImag) ivec.push_back(a.imag(), &inp_a.real());
    if (inp_b.hasReal() and bHasImag) ivec.push_back(b.imag(), &inp_b.real());
    if (inp_a.hasImag() and aHasReal) ivec.push_back(a_conj*a.real(), &inp_a.imag());
    if (inp_b.hasImag() and bHasReal) ivec.push_back(b_conj*b.real(), &inp_b.imag());

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

/*
void add(double prec, Orbital &out, ComplexVector &coefs, OrbitalVector &inp, bool union_grid) {
    NEEDS_TESTING;
    // set output spin from first contributing input
    for (int i = 0; i < inp.size(); i++) {
        if (abs(coefs[i]) < mrcpp::MachineZero) continue;
        if (inp[i]->getOccupancy() == 0) continue;
        out.setSpin(inp[i]->getSpin());
        break;
    }
    // all contributing input spins must be equal
    for (int i = 0; i < inp.size(); i++) {
        if (abs(coefs[i]) < mrcpp::MachineZero) continue;
        if (inp[i]->getOccupancy() == 0) continue;
        if (out.getSpin() != inp[i]->getSpin()) MSG_FATAL("Mixing spins");
    }

    add(prec, (QMFunction &) out, coefs, (QMFunctionVector &) inp, union_grid);
}
*/

/** out = inp_a * inp_b
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves. */
Orbital multiply(Orbital inp_a, Orbital inp_b, double prec) {
    int occ = compare_occ(inp_a, inp_b);
    int spin = compare_spin(inp_a, inp_b);
    Orbital out(occ, spin);

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

int compare_occ(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.occ() == orb_b.occ()) comp = orb_a.occ();
    return comp;
}

int compare_spin(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.spin() == orb_b.spin()) comp = orb_a.spin();
    return comp;
}

///** Normalize all orbitals in vector */
//void normalize(OrbitalVector &out) {
//    for (int i = 0; i < out.size(); i++) {
//        Orbital &out_i = *out[i];
//        qmfunction::normalize(out_i);
//    }
//}

} //namespace orbital


namespace density {

/****************************************
 * Density related standalone functions *
 ****************************************/

} //namespace density

} //namespace mrchem
