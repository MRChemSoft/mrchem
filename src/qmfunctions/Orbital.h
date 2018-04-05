#pragma once

#include "QMFunction.h"

/** @class Orbital
 *
 * @brief General complex-valued function with spin
 *
 * Inherits the general features of a complex function from QMFunction which
 * means separate MW function representations for the real and imaginary parts.
 * Note that there are several options for copying/assignment: the proper copy
 * constructor and assignment operator are *shallow* copies, which means that
 * they simply copy the *re and *im pointers (no transfer of ownership).
 * Additionaly, there is a deepCopy() which returns a *full* copy of the orbital,
 * and a paramCopy() which returns an empty orbital with the same rank_id/spin/occ.
 *
 * NOTE: since the standard copies are shallow copies and several orbitals can
 * point to the same MW functions, it is YOUR responibility to keep track of the
 * ownership and free the FunctionTree pointers before the last orbital goes out
 * of scope.
 */

namespace mrchem {

/* POD struct for orbital meta data. Used for simple MPI communication. */
struct OrbitalMeta {
    int rank_id;
    int spin;
    int occ;
    int nChunksReal;
    int nChunksImag;
    bool conjugate;
    double error;
};

class Orbital final : public QMFunction {
public:
    Orbital();
    Orbital(int spin, int occ = -1, int rank = -1);

    Orbital(const Orbital &orb);
    Orbital &operator=(const Orbital &orb);
    Orbital paramCopy() const;
    Orbital deepCopy();
    Orbital dagger() const;

    void setError(double error) { this->meta.error = error; }
    void setSpin(int spin) { this->meta.spin = spin; }
    void setOcc(int occ) { this->meta.occ = occ; }

    OrbitalMeta &getMetaData();
    int occ() const { return this->meta.occ; }
    int spin() const { return this->meta.spin; }
    int rankID() const { return this->meta.rank_id; }
    bool conjugate() const { return this->meta.conjugate; }
    double error() const { return this->meta.error; }
    double norm() const;
    double squaredNorm() const;

    void add(ComplexDouble c, Orbital inp, double prec = -1.0);
    void multiply(Orbital inp, double prec = -1.0);
    void rescale(ComplexDouble c);

    void normalize() { rescale(1.0/this->norm()); }
    void orthogonalize(Orbital inp);
    void orthogonalize(OrbitalVector inp_vec);

    void saveOrbital(const std::string &file);
    void loadOrbital(const std::string &file);

    char printSpin() const;
    friend std::ostream& operator<<(std::ostream &o, Orbital orb) { return orb.print(o); }

protected:
    OrbitalMeta meta;

    std::ostream& print(std::ostream &o) const;
};

} //namespace mrchem
