#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/mathutils.h"
#include "utils/RRMaximizer.h"

#include "qmfunctions.h"
#include "Orbital.h"
#include "Density.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace orbital {

/****************************************
 * Orbital related standalone functions *
 ****************************************/

/** Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Complicated by the fact that both bra and ket can be interpreted
 *  as complex conjugate versions of themselves. Notice that the <bra|
 *  position is already complex conjugated.
 *
 */
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

    return multiply(coefs, orbs, prec);
}

/** out_i = a*(inp_a)_i + b*(inp_b)_i
 *
 *  Component-wise addition of orbitals.
 *
 */
OrbitalVector add(ComplexDouble a, OrbitalVector &inp_a,
                  ComplexDouble b, OrbitalVector &inp_b,
                  double prec) {
    if (inp_a.size() != inp_b.size()) MSG_ERROR("Size mismatch");

    OrbitalVector out;
    for (int i = 0; i < inp_a.size(); i++) {
        Orbital out_i = add(a, inp_a[i], b, inp_b[i]);
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

/** out = c_0*inp_0 + c_1*inp_1 + ...
  *
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves.
  *
  */
Orbital multiply(const ComplexVector &c, OrbitalVector &inp, double prec) {
    if (c.size() != inp.size()) MSG_ERROR("Size mismatch");
    double thrs = mrcpp::MachineZero;

    Orbital out;
    // set output spin from first contributing input
    for (int i = 0; i < inp.size(); i++) {
        if (std::abs(c[i]) < thrs) continue;
        out = inp[i].paramCopy();
        break;
    }
    // all contributing input spins must be equal
    for (int i = 0; i < inp.size(); i++) {
        if (std::abs(c[i]) < thrs) continue;
        if (out.spin() != inp[i].spin()) {
            // empty orbitals with wrong spin can occur
            if (inp[i].hasReal()) MSG_FATAL("Mixing spins");
            if (inp[i].hasImag()) MSG_FATAL("Mixing spins");
        }
    }

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    for (int i = 0; i < inp.size(); i++) {
        bool cHasReal = (std::abs(c[i].real()) > thrs);
        bool cHasImag = (std::abs(c[i].imag()) > thrs);

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
OrbitalVector multiply(const ComplexMatrix &U, OrbitalVector &inp, double prec) {
    if (mpi::orb_size > 1) NOT_IMPLEMENTED_ABORT;

    OrbitalVector out;
    for (int i = 0; i < U.rows(); i++) {
        const ComplexVector &c = U.row(i);
        Orbital out_i = multiply(c, inp, prec);
        out.push_back(out_i);
    }
    return out;
}

/** Deep copy
 *
 * New orbitals are constructed as deep copies of the input set.
 *
 */
OrbitalVector deep_copy(OrbitalVector &inp) {
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = inp[i].deepCopy();
        out.push_back(out_i);
    }
    return out;
}

/** Parameter copy
 *
 * New orbitals are constructed as parameter copies of the input set.
 *
 */
OrbitalVector param_copy(const OrbitalVector &inp) {
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = inp[i].paramCopy();
        out.push_back(out_i);
    }
    return out;
}

/** Adjoin two vectors
 *
 * The orbitals of the input vector are appended to
 * (*this) vector, the ownership is transferred. Leaves
 * the input vector empty.
 *
 */
OrbitalVector adjoin(OrbitalVector &inp_a, OrbitalVector &inp_b) {
    OrbitalVector out;
    for (int i = 0; i < inp_a.size(); i++) out.push_back(inp_a[i]);
    for (int i = 0; i < inp_b.size(); i++) out.push_back(inp_b[i]);
    inp_a.clear();
    inp_b.clear();
    return out;
}

/** Disjoin vector in two parts
 *
 * All orbitals of a particular spin is collected in a new vector
 * and returned. These orbitals are removed from (*this) vector,
 * and the ownership is transferred.
 *
 */
OrbitalVector disjoin(OrbitalVector &inp, int spin) {
    OrbitalVector out;
    OrbitalVector tmp;
    for (int i = 0; i < inp.size(); i++) {
        Orbital &orb_i = inp[i];
        if (orb_i.spin() == spin) {
            out.push_back(orb_i);
        } else {
            tmp.push_back(orb_i);
        }
    }
    inp.clear();
    inp = tmp;
    return out;
}

/** Frees each orbital in the vector
 *
 * Leaves an empty vector. Orbitals are freed.
 *
 */
void free(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) vec[i].free();
    vec.clear();
}

/** Normalize all orbitals in the set */
void normalize(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i].normalize();
    }
}

/** Gram-Schmidt orthogonalize orbitals within the set */
void orthogonalize(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) {
        Orbital &orb_i = vec[i];
        for (int j = 0; j < i; j++) {
            Orbital &orb_j = vec[j];
            orb_i.orthogonalize(orb_j);
        }
    }
}

/** Orthogonalize the out orbital against all orbitals in inp */
void orthogonalize(OrbitalVector &vec, OrbitalVector &inp) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i].orthogonalize(inp);
    }
}

ComplexMatrix calc_overlap_matrix(OrbitalVector &braket) {
    return orbital::calc_overlap_matrix(braket, braket);
}

ComplexMatrix calc_overlap_matrix(OrbitalVector &bra, OrbitalVector &ket) {
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix S(Ni, Nj);
    S.setZero();

    for (int i = 0; i < Ni; i++) {
        for (int j = 0; j < Nj; j++) {
            S(i,j) = orbital::dot(bra[i], ket[j]);
        }
    }
    return S;
}

/** @brief Compute Löwdin orthonormalization matrix
 *
 * @param Phi: orbitals to orthonomalize
 *
 * Computes the inverse square root of the orbital overlap matrix S^(-1/2)
 */
ComplexMatrix calc_lowdin_matrix(OrbitalVector &Phi) {
    Timer timer;
    printout(1, "Calculating Löwdin orthonormalization matrix      ");

    ComplexMatrix S_tilde = orbital::calc_overlap_matrix(Phi);
    ComplexMatrix S_m12 = mathutils::hermitian_matrix_pow(S_tilde, -1.0/2.0);

    timer.stop();
    println(1, timer.getWallTime());
    return S_m12;
}

/** @brief Minimize the spatial extension of orbitals, by orbital rotation
 *
 * @param Phi: orbitals to localize
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 * which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix localize(double prec, OrbitalVector &Phi) {
    Printer::printHeader(0, "Localizing orbitals");
    Timer timer;

    ComplexMatrix U;
    int n_it = 0;
    if (Phi.size() > 1) {
        Timer rr_t;
        RRMaximizer rr(prec, Phi);
        n_it = rr.maximize();
        rr_t.stop();
        Printer::printDouble(0, "Computing Foster-Boys matrix", rr_t.getWallTime(), 5);

        if (n_it > 0) {
            println(0, " Converged after iteration   " << std::setw(30) << n_it);
            U = rr.getTotalU().transpose().cast<ComplexDouble>();
        } else {
            println(0, " Foster-Boys localization did not converge!");
        }
    } else {
        println(0, " Cannot localize less than two orbitals");
    }

    if (n_it <= 0) {
        Timer orth_t;
        U = orbital::calc_lowdin_matrix(Phi);
        orth_t.stop();
        Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);
    }

    Timer rot_t;
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return U;
}

/** @brief Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * @param Phi: orbitals to rotate
 * @param F: Fock matrix to diagonalize
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place and Fock matrix is diagonalized in place.
 * The transformation matrix is returned.
 */
ComplexMatrix diagonalize(double prec, OrbitalVector &Phi, ComplexMatrix &F) {
    Printer::printHeader(0, "Digonalizing Fock matrix");
    Timer timer;

    Timer orth_t;
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    F = S_m12.transpose()*F*S_m12;
    orth_t.stop();
    Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);

    Timer diag_t;
    ComplexMatrix U = ComplexMatrix::Zero(F.rows(), F.cols());
    int np = orbital::size_paired(Phi);
    int na = orbital::size_alpha(Phi);
    int nb = orbital::size_beta(Phi);
    if (np > 0) mathutils::diagonalize_block(F, U, 0,       np);
    if (na > 0) mathutils::diagonalize_block(F, U, np,      na);
    if (nb > 0) mathutils::diagonalize_block(F, U, np + na, nb);
    U = U * S_m12;
    diag_t.stop();
    Printer::printDouble(0, "Diagonalizing matrix", diag_t.getWallTime(), 5);

    Timer rot_t;
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return U;
}

/** @brief Perform the Löwdin orthonormalization
 *
 * @param Phi: orbitals to orthonormalize
 *
 * Orthonormalizes the orbitals by multiplication of the Löwdin matrix S^(-1/2).
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix orthonormalize(double prec, OrbitalVector &Phi) {
    ComplexMatrix U = orbital::calc_lowdin_matrix(Phi);
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    return U;
}

/** Returns the number of occupied orbitals */
int size_occupied(const OrbitalVector &vec) {
    int nOcc = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() > 0) nOcc++;
    }
    return nOcc;
}

/** Returns the number of empty orbitals */
int size_empty(const OrbitalVector &vec) {
    int nEmpty = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 0) nEmpty++;
    }
    return nEmpty;
}

/** Returns the number of singly occupied orbitals */
int size_singly(const OrbitalVector &vec) {
    int nSingly = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 1) nSingly++;
    }
    return nSingly;
}

/** Returns the number of doubly occupied orbitals */
int size_doubly(const OrbitalVector &vec) {
    int nDoubly = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 1) nDoubly++;
    }
    return nDoubly;
}

/** Returns the number of paired orbitals */
int size_paired(const OrbitalVector &vec) {
    int nPaired = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Paired) nPaired++;
    }
    return nPaired;
}

/** Returns the number of alpha orbitals */
int size_alpha(const OrbitalVector &vec) {
    int nAlpha = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Alpha) nAlpha++;
    }
    return nAlpha;
}

/** Returns the number of beta orbitals */
int size_beta(const OrbitalVector &vec) {
    int nBeta = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Beta) nBeta++;
    }
    return nBeta;
}

/** Returns the spin multiplicity of the vector */
int get_multiplicity(const OrbitalVector &vec) {
    int nAlpha = get_electron_number(vec, SPIN::Alpha);
    int nBeta = get_electron_number(vec, SPIN::Beta);
    int S = std::abs(nAlpha - nBeta);
    return S + 1;
}

/** Returns the number of electrons with the given spin
 *
 * Paired spin (default input) returns the total number of electrons.
 *
 */
int get_electron_number(const OrbitalVector &vec, int spin) {
    int nElectrons = 0;
    for (int i = 0; i < vec.size(); i++) {
        const Orbital &orb = vec[i];
        if (spin == SPIN::Paired) {
            nElectrons += orb.occ();
        } else if (spin == SPIN::Alpha) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Alpha) {
                nElectrons += 1;
            }
        } else if (spin == SPIN::Beta) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Beta) {
                nElectrons += 1;
            }
        } else {
            MSG_ERROR("Invalid spin argument");
        }
    }
    return nElectrons;
}

/** Returns a vector containing the orbital errors */
DoubleVector get_errors(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector errors = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        errors(i) = vec[i].error();
    }
    return errors;
}

/** Assign errors to each orbital.
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_errors(OrbitalVector &vec, const DoubleVector &errors) {
    if (vec.size() != errors.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setError(errors(i));
    }
}

/** Returns a vector containing the orbital spins */
IntVector get_spins(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    IntVector spins = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        spins(i) = vec[i].spin();
    }
    return spins;
}

/** Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_spins(OrbitalVector &vec, const IntVector &spins) {
    if (vec.size() != spins.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setSpin(spins(i));
    }
}

/** Returns a vector containing the orbital occupancies */
IntVector get_occupancies(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    IntVector occ = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        occ(i) = vec[i].occ();
    }
    return occ;
}

/** Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_occupancies(OrbitalVector &vec, const IntVector &occ) {
    if (vec.size() != occ.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setOcc(occ(i));
    }
}

/** Returns a vector containing the orbital square norms */
DoubleVector get_squared_norms(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(vec[i])) norms(i) = vec[i].squaredNorm();
    }
    return norms;
}

/** Returns a vector containing the orbital norms */
DoubleVector get_norms(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(vec[i])) norms(i) = vec[i].norm();
    }
    return norms;
}

} //namespace orbital


namespace density {

/****************************************
 * Density related standalone functions *
 ****************************************/
void calc_density(Density &rho, Orbital phi, double prec) {
    if (rho.hasTotal()) MSG_ERROR("Density not empty");
    if (rho.hasSpin())  MSG_ERROR("Density not empty");
    if (rho.hasAlpha()) MSG_ERROR("Density not empty");
    if (rho.hasBeta())  MSG_ERROR("Density not empty");

    double occ = 1.0;
    if (not rho.isSpinDensity()) occ = (double) phi.occ();

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*real_2, phi.real());
        mrcpp::multiply(prec, *real_2, occ, phi.real(), phi.real(), 1);
        sum_vec.push_back(real_2);
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*imag_2, phi.imag());
        mrcpp::multiply(prec, *imag_2, occ, phi.imag(), phi.imag(), 1);
        sum_vec.push_back(imag_2);
    }

    if (rho.isSpinDensity()) {
        if (phi.spin() == SPIN::Paired) {
            rho.alloc(DENSITY::Alpha);
            mrcpp::copy_grid(rho.alpha(), sum_vec);
            mrcpp::add(-1.0, rho.alpha(), sum_vec, 0);

            rho.alloc(DENSITY::Beta);
            mrcpp::copy_grid(rho.beta(), sum_vec);
            mrcpp::add(-1.0, rho.beta(), sum_vec, 0);
        }
        if (phi.spin() == SPIN::Alpha) {
            rho.alloc(DENSITY::Alpha);
            mrcpp::copy_grid(rho.alpha(), sum_vec);
            mrcpp::add(-1.0, rho.alpha(), sum_vec, 0);

            rho.alloc(DENSITY::Beta);
            mrcpp::copy_grid(rho.beta(), sum_vec);
            rho.beta().setZero();
        }
        if (phi.spin() == SPIN::Beta) {
            rho.alloc(DENSITY::Alpha);
            mrcpp::copy_grid(rho.alpha(), sum_vec);
            rho.alpha().setZero();

            rho.alloc(DENSITY::Beta);
            mrcpp::copy_grid(rho.beta(), sum_vec);
            mrcpp::add(-1.0, rho.beta(), sum_vec, 0);
        }
        FunctionTreeVector<3> tot_vec;
        tot_vec.push_back(1.0, &rho.alpha());
        tot_vec.push_back(1.0, &rho.beta());
        rho.alloc(DENSITY::Total);
        mrcpp::copy_grid(rho.total(), tot_vec);
        mrcpp::add(-1.0, rho.total(), tot_vec, 0);

        FunctionTreeVector<3> spin_vec;
        spin_vec.push_back( 1.0, &rho.alpha());
        spin_vec.push_back(-1.0, &rho.beta());
        rho.alloc(DENSITY::Spin);
        mrcpp::copy_grid(rho.spin(), spin_vec);
        mrcpp::add(-1.0, rho.spin(), spin_vec, 0);
    } else {
        rho.alloc(DENSITY::Total);
        mrcpp::copy_grid(rho.total(), sum_vec);
        mrcpp::add(-1.0, rho.total(), sum_vec, 0);
    }
    sum_vec.clear(true);
}

void calc_density(Density &rho, OrbitalVector &Phi, double prec) {
    if (rho.hasTotal()) MSG_ERROR("Density not empty");
    if (rho.hasSpin())  MSG_ERROR("Density not empty");
    if (rho.hasAlpha()) MSG_ERROR("Density not empty");
    if (rho.hasBeta())  MSG_ERROR("Density not empty");

    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i

    FunctionTreeVector<3> total_vec, spin_vec, alpha_vec, beta_vec;
    Density* rho_tmp1 = 0;
    Density* rho_tmp2 = 0;
    Density* rho_i = 0;
    if (mpi::orb_size > 1) {
        //only master does the summation
        for (int i_orb = 0; i_orb < Phi.size(); i_orb++) {
            rho_i = new Density(rho);
            if (i_orb%mpi::orb_size == mpi::orb_rank) {
                Orbital &phi_i = Phi[i_orb];
                calc_density(*rho_i, phi_i, mult_prec);
                if (mpi::orb_rank != 0) mpi::send_density(*rho_i, 0, i_orb);
            }
            //add on the fly
            if (mpi::orb_rank == 0 and i_orb == 0){
                //first iteration does not sum only receive in tmp1
                if (i_orb%mpi::orb_size != mpi::orb_rank) {
                    mpi::recv_density(*rho_tmp1, i_orb%mpi::orb_size, i_orb);
                } else {
                    rho_tmp1 = rho_i;
                    rho_i = rho_tmp2;
                }
            } else if (mpi::orb_rank == 0) {
                if (i_orb%mpi::orb_size != mpi::orb_rank) {
                    mpi::recv_density(*rho_i, i_orb%mpi::orb_size, i_orb);
                }
                //exchange pointers. old result is in tmp1: tmp1 = rho_i+tmp2
                rho_tmp2 = rho_tmp1;
                rho_tmp1 = new Density(rho);
                if (i_orb == Phi.size() - 1) {
                    //last iteration, put result into rho
                    rho_tmp1 = &rho;
                }
                if (rho_i->hasTotal()) {
                    if (not rho_tmp1->hasTotal()) rho_tmp1->alloc(DENSITY::Total);
                    total_vec.push_back(&rho_tmp2->total());
                    total_vec.push_back(&rho_i->total());
                    mrcpp::copy_grid(rho_tmp1->total(), total_vec);
                    mrcpp::add(-1.0, rho_tmp1->total(), total_vec, 0);
                    total_vec.clear(true);
                }
                if (rho_i->hasSpin()) {
                    if (not rho_tmp1->hasSpin()) rho_tmp1->alloc(DENSITY::Spin);
                    spin_vec.push_back(&rho_tmp2->spin());
                    spin_vec.push_back(&rho_i->spin());
                    mrcpp::copy_grid(rho_tmp1->spin(), spin_vec);
                    mrcpp::add(-1.0, rho_tmp1->spin(), spin_vec, 0);
                    spin_vec.clear(true);
                }
                if (rho_i->hasAlpha()) {
                    if (not rho_tmp1->hasAlpha()) rho_tmp1->alloc(DENSITY::Alpha);
                    alpha_vec.push_back(&rho_tmp2->alpha());
                    alpha_vec.push_back(&rho_i->alpha());
                    mrcpp::copy_grid(rho_tmp1->alpha(), alpha_vec);
                    mrcpp::add(-1.0, rho_tmp1->alpha(), alpha_vec, 0);
                    alpha_vec.clear(true);
                }
                if (rho_i->hasBeta()) {
                    if (not rho_tmp1->hasBeta()) rho_tmp1->alloc(DENSITY::Beta);
                    beta_vec.push_back(&rho_tmp2->beta());
                    beta_vec.push_back(&rho_i->beta());
                    mrcpp::copy_grid(rho_tmp1->beta(), beta_vec);
                    mrcpp::add(-1.0, rho_tmp1->beta(), beta_vec, 0);
                    beta_vec.clear(true);
                }
            }
        }
    } else {
        //Serial processing
        for (int i = 0; i < Phi.size(); i++) {
            Density rho_i(rho);
            calc_density(rho_i, Phi[i], mult_prec);
            if (rho_i.hasTotal()) total_vec.push_back(&rho_i.total());
            if (rho_i.hasSpin())   spin_vec.push_back(&rho_i.spin());
            if (rho_i.hasAlpha()) alpha_vec.push_back(&rho_i.alpha());
            if (rho_i.hasBeta())   beta_vec.push_back(&rho_i.beta());
            rho_i.clear(); // not free
        }
    }

    if (mpi::orb_rank == 0 and not rho.isShared()) {
        if (total_vec.size() > 5) {
            rho.alloc(DENSITY::Total);
            mrcpp::add(add_prec, rho.total(), total_vec);
        } else if (total_vec.size() > 0) {
            rho.alloc(DENSITY::Total);
            mrcpp::copy_grid(rho.total(), total_vec);
            mrcpp::add(-1.0, rho.total(), total_vec, 0);
        }
        if (spin_vec.size() > 5) {
            rho.alloc(DENSITY::Spin);
            mrcpp::add(add_prec, rho.spin(), spin_vec);
        } else if (spin_vec.size() > 0) {
            rho.alloc(DENSITY::Spin);
            mrcpp::copy_grid(rho.spin(), spin_vec);
            mrcpp::add(-1.0, rho.spin(), spin_vec, 0);
        }
        if (alpha_vec.size() > 5) {
            rho.alloc(DENSITY::Alpha);
            mrcpp::add(add_prec, rho.alpha(), alpha_vec);
        } else if (alpha_vec.size() > 0) {
            rho.alloc(DENSITY::Alpha);
            mrcpp::copy_grid(rho.alpha(), alpha_vec);
            mrcpp::add(-1.0, rho.alpha(), alpha_vec, 0);
        }
        if (beta_vec.size() > 5) {
            rho.alloc(DENSITY::Beta);
            mrcpp::add(add_prec, rho.beta(), beta_vec);
        } else if (beta_vec.size() > 0) {
            rho.alloc(DENSITY::Beta);
            mrcpp::copy_grid(rho.beta(), beta_vec);
            mrcpp::add(-1.0, rho.beta(), beta_vec, 0);
        }
        total_vec.clear(true);
        spin_vec.clear(true);
        alpha_vec.clear(true);
        beta_vec.clear(true);
    }

    if (mpi::orb_size > 1) {
        //we always broadcast density
        //If the density is shared, only metadata will be sent/received
        if (mpi::orb_rank == 0) {
            for (int i_MPI = 1; i_MPI < mpi::orb_size; i_MPI++) {
                mpi::send_density(rho, i_MPI, 54);
            }
        } else {
            //do nothing, only receive Density
            mpi::recv_density(rho, 0, 54);
        }
    }
}

} //namespace density

} //namespace mrchem

/** Collective printing of all orbitals in the vector */
std::ostream& operator<<(std::ostream &o, const mrchem::OrbitalVector &vec) {
    Printer::setScientific();
    o << "*OrbitalVector: ";
    o << std::setw(4) << vec.size()          << " orbitals  ";
    o << std::setw(4) << mrchem::orbital::size_occupied(vec)  << " occupied  ";
    o << std::setw(4) << mrchem::orbital::get_electron_number(vec) << " electrons " << std::endl;
    o << "------------------------------------------------------------\n";
    o << "   n  RankID           Norm          Spin Occ      Error    \n";
    o << "------------------------------------------------------------\n";
    for (int i = 0; i < vec.size(); i++) {
        o << std::setw(4) << i << vec[i] << std::endl;
    }
    o << "------------------------------------------------------------\n";
    return o;
}
