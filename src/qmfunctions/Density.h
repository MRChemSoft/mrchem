#pragma once

#include "QMFunction.h"
#include "density_utils.h"

/** @class Density
 *
 * @brief General complex-valued function to handle densities (incl trans. densities)
 *
 * Inherits the general features of a complex function from QMFunction which
 * means separate MW function representations for the real and imaginary parts.
 * Note that there are several options for copying/assignment: the proper copy
 * constructor and assignment operator are *shallow* copies, which means that
 * they simply copy the *re and *im pointers (no transfer of ownership).
 * Additionaly, there is a deepCopy() which returns a *full* copy of the density,
 * and a paramCopy() which returns an empty density with the same rank_id/spin.
 *
 * NOTE: since the standard copies are shallow copies and several densitys can
 * point to the same MW functions, it is YOUR responibility to keep track of the
 * ownership and free the FunctionTree pointers before the last density goes out
 * of scope.
 */

namespace mrchem {

/* POD struct for density meta data. Used for simple MPI communication. */
struct DensityMeta {
    int rank_id;
    int spin;
    int nChunksReal;
    int nChunksImag;
    bool conjugate;
    double error;
};

class Density final : public QMFunction {
public:
    Density();
    Density(int spin, int rank = -1);

    Density(const Density &dens);
    Density &operator=(const Density &dens);
    Density paramCopy() const;
    Density deepCopy();
    Density dagger() const;
    
    void setError(double error) { this->meta.error = error; }
    void setRankId(int rank) { this->meta.rank_id = rank; }
    void setSpin(int spin) { this->meta.spin = spin; }

    void allocReal();
    void allocImag();
    
    DensityMeta &getMetaData();
    int spin() const { return this->meta.spin; }
    int rankID() const { return this->meta.rank_id; }
    bool conjugate() const { return this->meta.conjugate; }
    double error() const { return this->meta.error; }
    double norm() const;
    double squaredNorm() const;

    void multiply(Density inp, double prec = -1.0);
    
    void saveDensity(const std::string &file);
    void loadDensity(const std::string &file);

    char printSpin() const;
    friend std::ostream& operator<<(std::ostream &o, Density dens) { return dens.print(o); }

protected:
    DensityMeta meta;

    std::ostream& print(std::ostream &o) const;
};

} //namespace mrchem
