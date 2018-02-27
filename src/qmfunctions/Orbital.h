#pragma once

#include "qmfunctions.h"

namespace mrchem {

/* POD struct for orbital meta data. Used for simple MPI communication. */
struct OrbitalMeta {
    int occ;
    int spin;
    int nChunksReal;
    int nChunksImag;
    bool conjugate;
    double error;
};

class Orbital {
public:
    Orbital(int occ = 0, int s = 0);
    ~Orbital() { }

    Orbital(const Orbital &orb);
    Orbital &operator=(const Orbital &orb);
    Orbital deepCopy();
    Orbital dagger() const;

    void alloc(int type = NUMBER::Total);
    void clear(int type = NUMBER::Total);
    void free(int type = NUMBER::Total);

    bool hasReal() const { return (this->re == 0) ? false : true; }
    bool hasImag() const { return (this->im == 0) ? false : true; }

    mrcpp::FunctionTree<3> &real() { return *this->re; }
    mrcpp::FunctionTree<3> &imag() { return *this->im; }

    const mrcpp::FunctionTree<3> &real() const { return *this->re; }
    const mrcpp::FunctionTree<3> &imag() const { return *this->im; }

    void setReal(mrcpp::FunctionTree<3> *real) { this->re = real; }
    void setImag(mrcpp::FunctionTree<3> *imag) { this->im = imag; }

    void setError(double err) { this->meta.error = err; }

    OrbitalMeta &getMetaData();
    int occ() const { return this->meta.occ; }
    int spin() const { return this->meta.spin; }
    bool conjugate() const { return this->meta.conjugate; }
    double error() const { return this->meta.error; }
    double norm() const;
    double squaredNorm() const;

    void add(ComplexDouble c, Orbital inp, double prec = -1.0);
    void multiply(Orbital inp, double prec = -1.0);
    void rescale(ComplexDouble c);

    void normalize() { rescale(1.0/this->norm()); }
    void orthogonalize(Orbital inp);

    friend std::ostream& operator<<(std::ostream &o, Orbital orb) { return orb.print(o); }

protected:
    OrbitalMeta meta;
    mrcpp::FunctionTree<3> *re;     ///* Real part of function
    mrcpp::FunctionTree<3> *im;     ///* Imaginary part of function

    char printSpin() const;
    std::ostream& print(std::ostream &o) const;
};

} //namespace mrchem
