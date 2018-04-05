#include <fstream>

#include "MRCPP/Printer"

#include "Orbital.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief Default constructor
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 */
Orbital::Orbital()
        : QMFunction(0, 0),
          meta({-1, 0, 0, 0, 0, false, 0.0}) {
}

/** @brief Constructor
 *
 * @param spin: electron spin (SPIN::Alpha/Beta/Paired)
 * @param occ: occupancy
 * @param rank: MPI ownership (-1 means all MPI ranks)
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 */
Orbital::Orbital(int spin, int occ, int rank)
        : QMFunction(0, 0),
          meta({rank, spin, occ, 0, 0, false, 0.0}) {
    if (this->spin() < 0) INVALID_ARG_ABORT;
    if (this->occ() < 0) {
        if (this->spin() == SPIN::Paired) this->meta.occ = 2;
        if (this->spin() == SPIN::Alpha) this->meta.occ = 1;
        if (this->spin() == SPIN::Beta) this->meta.occ = 1;
    }
}

/** @brief Copy constructor
 *
 * @param orb: orbital to copy
 *
 * Shallow copy: meta data is copied along with the *re and *im pointers,
 * NO transfer of ownership.
 */
Orbital::Orbital(const Orbital &orb)
        : QMFunction(orb),
          meta(orb.meta) {
}

/** @brief Assignment operator
 *
 * @param orb: orbital to copy
 *
 * Shallow copy: meta data is copied along with the *re and *im pointers,
 * NO transfer of ownership.
 */
Orbital& Orbital::operator=(const Orbital &orb) {
    if (this != &orb) {
        this->re = orb.re;
        this->im = orb.im;
        this->meta = orb.meta;
    }
    return *this;
}

/** @brief Parameter copy
 *
 * Returns a new orbital with the same spin, occupancy and rank_id as *this orbital.
 */
Orbital Orbital::paramCopy() const {
    return Orbital(this->spin(), this->occ(), this->rankID());
}

/** @brief Deep copy
 *
 * Returns a new orbital which is a full blueprint copy of *this orbital. This is
 * achieved by building a new grid for the real and imaginary parts and adding
 * in place.
 */
Orbital Orbital::deepCopy() {
    Orbital out(*this); // Shallow copy (should copy all meta data)
    out.clear();        // Remove *re and *im pointers
    if (this->hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::FunctionTreeVector<3> vec;
        vec.push_back(&this->real());
        mrcpp::copy_grid(out.real(), vec);
        mrcpp::add(-1.0, out.real(), vec);
    }
    if (this->hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::FunctionTreeVector<3> vec;
        vec.push_back(&this->imag());
        mrcpp::copy_grid(out.imag(), vec);
        mrcpp::add(-1.0, out.imag(), vec);
    }
    return out;         // Return shallow copy
}

/** @brief Complex conjugation
 *
 * Returns a new orbital which is a shallow copy of *this orbital, with a flipped
 * conjugate parameter. Pointer ownership is not transferred, so *this and the output
 * orbital points to the same MW representations of the real and imaginary parts,
 * however, they interpret the imaginary part with opposite sign.
 */
Orbital Orbital::dagger() const {
    Orbital out(*this); // Shallow copy
    out.meta.conjugate = not this->meta.conjugate;
    return out;         // Return shallow copy
}

/** @brief Returns the orbital meta data
 *
 * Tree sizes (nChunks) are flushed before return.
 */
OrbitalMeta& Orbital::getMetaData() {
    this->meta.nChunksReal = 0;
    this->meta.nChunksImag = 0;
    if (this->hasReal()) this->meta.nChunksReal = real().getNChunksUsed();
    if (this->hasImag()) this->meta.nChunksImag = imag().getNChunksUsed();
    return this->meta;
}

/** @brief Returns the norm of the orbital */
double Orbital::norm() const {
    double norm = this->squaredNorm();
    if (norm > 0.0) norm = sqrt(norm);
    return norm;
}

/** @brief Returns the squared norm of the orbital */
double Orbital::squaredNorm() const {
    double sq_r = -1.0;
    double sq_i = -1.0;
    if (this->hasReal()) sq_r = this->real().getSquareNorm();
    if (this->hasImag()) sq_i = this->imag().getSquareNorm();

    double sq_norm = 0.0;
    if (sq_r < 0.0 and sq_i < 0.0) {
        sq_norm = -1.0;
    } else {
        if (sq_r >= 0.0) sq_norm += sq_r;
        if (sq_i >= 0.0) sq_norm += sq_i;
    }
    return sq_norm;
}

/** @brief In place multiply with scalar */
void Orbital::rescale(ComplexDouble c) {
    double thrs = mrcpp::MachineZero;
    bool cHasReal = (std::abs(c.real()) > thrs);
    bool cHasImag = (std::abs(c.imag()) > thrs);

    if (cHasReal and cHasImag) {
        Orbital tmp = orbital::add(c, *this, 0.0, *this);
        this->free();
        *this = tmp;
    }
    if (cHasReal and not cHasImag) {
        if (this->hasReal()) this->real().rescale(c.real());
        if (this->hasImag()) this->imag().rescale(c.real());
    }
    if (not cHasReal and not cHasImag) {
        if (this->hasReal()) this->real().setZero();
        if (this->hasImag()) this->imag().setZero();
    }
    if (not cHasReal and cHasImag) {
        double conj = 1.0;
        if (this->conjugate()) conj = -1.0;
        mrcpp::FunctionTree<3> *tmp_re = this->re;
        mrcpp::FunctionTree<3> *tmp_im = this->im;
        if (tmp_re != 0) tmp_re->rescale(c.imag());
        if (tmp_im != 0) tmp_im->rescale(-1.0*conj*c.imag());
        this->clear();
        this->setReal(tmp_im);
        this->setImag(tmp_re);
    }
}

/** @brief In place orthogonalize against inp */
void Orbital::orthogonalize(Orbital inp) {
    ComplexDouble overlap = orbital::dot(inp, *this);
    double sq_norm = inp.squaredNorm();
    if (std::abs(overlap) > mrcpp::MachineZero) {
        this->add(-1.0*(overlap/sq_norm), inp);
    }
}

/** @brief In place orthogonalize against all orbitals in inp vector */
void Orbital::orthogonalize(OrbitalVector inp_vec) {
    for (int i = 0; i < inp_vec.size(); i++) {
        this->orthogonalize(inp_vec[i]);
    }
}

/** @brief In place addition */
void Orbital::add(ComplexDouble c, Orbital inp, double prec) {
    // The following sets the spin/occupancy only locally, e.i. the
    // original orbital is NOT changed. The effect of this is to
    // disable spin/occupancy sanity checks in the following addition.
    // Since this is an in-place operation, we assume that the spin
    // and occupancy of the result is already defined and that it
    // should not change.
    inp.setSpin(this->spin());
    inp.setOcc(this->occ());
    Orbital tmp = orbital::add(1.0, *this, c, inp, prec);
    this->free();
    *this = tmp;
}

void Orbital::saveOrbital(const std::string &file) {
    //writing meta data
    std::stringstream metafile;
    metafile << file << ".meta";

    //this flushes tree sizes
    OrbitalMeta &my_meta = getMetaData();

    std::fstream f;
    f.open(metafile.str(), std::ios::out | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");
    f.write((char *) &my_meta, sizeof(OrbitalMeta));
    f.close();

    //writing real part
    if (hasReal()) {
        std::stringstream fname;
        fname << file << "_re";
        this->real().saveTree(fname.str());
    }

    //writing imaginary part
    if (hasImag()) {
        std::stringstream fname;
        fname << file << "_im";
        this->imag().saveTree(fname.str());
    }
}

void Orbital::loadOrbital(const std::string &file) {
    if (hasReal()) MSG_ERROR("Orbital not empty");
    if (hasImag()) MSG_ERROR("Orbital not empty");

    //reading meta data
    std::stringstream fmeta;
    fmeta << file << ".meta";

    //this flushes tree sizes
    OrbitalMeta &my_meta = getMetaData();

    std::fstream f;
    f.open(fmeta.str(), std::ios::in | std::ios::binary);
    if (f.is_open()) f.read((char *) &my_meta, sizeof(OrbitalMeta));
    f.close();

    //reading real part
    if (meta.nChunksReal > 0) {
        std::stringstream fname;
        fname << file << "_re";
        alloc(NUMBER::Real);
        this->real().loadTree(fname.str());
    }

    //reading imaginary part
    if (meta.nChunksImag > 0) {
        std::stringstream fname;
        fname << file << "_im";
        alloc(NUMBER::Imag);
        this->imag().loadTree(fname.str());
    }
}

/** @brief Returns a character representing the spin (a/b/p) */
char Orbital::printSpin() const {
    char sp = 'u';
    if (this->spin() == SPIN::Paired) sp = 'p';
    if (this->spin() == SPIN::Alpha) sp = 'a';
    if (this->spin() == SPIN::Beta) sp = 'b';
    return sp;
}

/** @brief Pretty output of orbital meta data */
std::ostream& Orbital::print(std::ostream &o) const {
    int oldprec = mrcpp::Printer::setPrecision(12);
    o << std::setw(6)  << this->rankID();
    o << std::setw(25) << this->norm();
    o << std::setw(5)  << this->printSpin();
    o << std::setw(4)  << this->occ();
    mrcpp::Printer::setPrecision(5);
    o << std::setw(15) << this->error();
    mrcpp::Printer::setPrecision(oldprec);
    return o;
}

} //namespace mrchem

