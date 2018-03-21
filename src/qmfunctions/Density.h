#pragma once

#include "qmfunctions.h"

namespace mrchem {

/** @class Density
 *
 *  @brief Placeholder for FunctionTrees that represent a density
 *
 *  A density contains FunctionTree pointers to each of the components
 *  (total, spin, alpha and beta). If the density is NOT a spin density,
 *  only the total component is used, otherwise ALL components are computed.
 *
 *  Density related stand-alone functions are located in:
 *  @file qmfunctions.h (for calculation of densities from orbitals)
 *  @file parallel.h (for MPI send/recv functionality)
 *
 */

/* POD struct for density meta data. Used for simple MPI communication. */
struct DensityMeta {
    bool is_spin;
    int nChunksTotal;
    int nChunksSpin;
    int nChunksAlpha;
    int nChunksBeta;
};

class Density final {
public:
    Density(bool spin, bool shared);
    Density(const Density &inp);
    Density& operator=(const Density &inp);
    ~Density();

    void alloc(int type = -1);
    void clear(int type = -1);
    void free(int type = -1);

    int getNNodes(int type = -1) const;
    DensityMeta &getMetaData();

    bool isSpinDensity() const { return this->meta.is_spin; }
    bool isShared() const { return (this->sh_mem == 0) ? false : true; }

    bool hasTotal() const { return (this->dens_t == 0) ? false : true; }
    bool hasSpin()  const { return (this->dens_s == 0) ? false : true; }
    bool hasAlpha() const { return (this->dens_a == 0) ? false : true; }
    bool hasBeta()  const { return (this->dens_b == 0) ? false : true; }

    void setDensity(mrcpp::FunctionTree<3> *rho, int type);

    mrcpp::FunctionTree<3> &total() { return *this->dens_t; }
    mrcpp::FunctionTree<3> &spin()  { return *this->dens_s; }
    mrcpp::FunctionTree<3> &alpha() { return *this->dens_a; }
    mrcpp::FunctionTree<3> &beta()  { return *this->dens_b; }

    const mrcpp::FunctionTree<3> &total() const { return *this->dens_t; }
    const mrcpp::FunctionTree<3> &spin()  const { return *this->dens_s; }
    const mrcpp::FunctionTree<3> &alpha() const { return *this->dens_a; }
    const mrcpp::FunctionTree<3> &beta()  const { return *this->dens_b; }

protected:
    DensityMeta meta;
    mrcpp::SharedMemory *sh_mem;
    mrcpp::FunctionTree<3> *dens_t;
    mrcpp::FunctionTree<3> *dens_s;
    mrcpp::FunctionTree<3> *dens_a;
    mrcpp::FunctionTree<3> *dens_b;
};

} //namespace mrchem
