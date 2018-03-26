#include "MRCPP/Printer"

#include "Orbital.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

Orbital::Orbital()
        : QMFunction(0, 0),
          meta({-1, 0, 0, 0, 0, false, 0.0}) {
}

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

Orbital::Orbital(const Orbital &orb)
        : QMFunction(orb),
          meta(orb.meta) {
}

Orbital& Orbital::operator=(const Orbital &orb) {
    if (this != &orb) {
        this->re = orb.re;
        this->im = orb.im;
        this->meta = orb.meta;
    }
    return *this;
}

Orbital Orbital::paramCopy() const {
    return Orbital(this->spin(), this->occ(), this->rankID());
}

Orbital Orbital::deepCopy() {
    Orbital out(*this); // Shallow copy
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

Orbital Orbital::dagger() const {
    Orbital out(*this); // Shallow copy
    out.meta.conjugate = not this->meta.conjugate;
    return out;         // Return shallow copy
}

/** Tree sizes (nChunks) are flushed before return. */
OrbitalMeta& Orbital::getMetaData() {
    this->meta.nChunksReal = 0;
    this->meta.nChunksImag = 0;
    if (this->hasReal()) this->meta.nChunksReal = real().getNChunksUsed();
    if (this->hasImag()) this->meta.nChunksImag = imag().getNChunksUsed();
    return this->meta;
}

double Orbital::norm() const {
    double norm = this->squaredNorm();
    if (norm > 0.0) norm = sqrt(norm);
    return norm;
}

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

/** In place multiply with scalar */
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

/** In place orthogonalize against inp */
void Orbital::orthogonalize(Orbital inp) {
    ComplexDouble overlap = orbital::dot(inp, *this);
    double sq_norm = inp.squaredNorm();
    if (std::abs(overlap) > mrcpp::MachineZero) {
        this->add(-1.0*(overlap/sq_norm), inp);
    }
}

/** In place orthogonalize against all orbitals in inp */
void Orbital::orthogonalize(OrbitalVector inp_vec) {
    for (int i = 0; i < inp_vec.size(); i++) {
        this->orthogonalize(inp_vec[i]);
    }
}

/** In place addition */
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

char Orbital::printSpin() const {
    char sp = 'u';
    if (this->spin() == SPIN::Paired) sp = 'p';
    if (this->spin() == SPIN::Alpha) sp = 'a';
    if (this->spin() == SPIN::Beta) sp = 'b';
    return sp;
}

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

