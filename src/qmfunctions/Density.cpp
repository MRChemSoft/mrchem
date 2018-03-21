#include "MRCPP/Printer"

#include "Density.h"
#include "parallel.h"

using mrcpp::FunctionTree;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param spin: compute spin densities
 * @param shared: use shared memory MPI
 *
 * If spin is activated all density components will be computed, otherwise only
 * total density is computed. If shared memory is activated all density components
 * will use the same shared memory block. Shared memory is only activated for MPI
 * runs with more than 1 process.
 */
Density::Density(bool spin, bool shared)
        : meta({spin, 0, 0, 0, 0}),
          sh_mem(0),
          dens_t(0),
          dens_s(0),
          dens_a(0),
          dens_b(0) {
    if (shared and mpi::share_size > 1) {
        //initiate up to 2000MB shared memory
        this->sh_mem = new mrcpp::SharedMemory(mpi::comm_share, 2000);
    }
}

/** @brief copy constructor
 *
 * @param inp: density to be copied
 *
 * Takes only the spin parameter from the input density, NOT shared memory info
 * and NOT the function data. All components are NULL.
 */
Density::Density(const Density &inp)
        : meta({inp.meta.is_spin, 0, 0, 0, 0}),
          sh_mem(0),
          dens_t(0),
          dens_s(0),
          dens_a(0),
          dens_b(0) {
}

/** @brief assignement operator
 *
 * @param inp: density to be copied
 *
 * Currently deactivated
 */
Density& Density::operator=(const Density &inp) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief destructor
 *
 * All density components should be cleared already at this point using
 * clear() or free().
 */
Density::~Density() {
    if (hasTotal()) MSG_ERROR("Density not properly deallocated");
    if (hasSpin())  MSG_ERROR("Density not properly deallocated");
    if (hasAlpha()) MSG_ERROR("Density not properly deallocated");
    if (hasBeta())  MSG_ERROR("Density not properly deallocated");
    if (isShared()) delete this->sh_mem;
}

/** @brief allocate density pointers
 *
 * @param type: density component
 *
 * This will allocate the FunctionTree pointer of the given density component.
 * The new FunctionTree will be uninitialized. All density components use the
 * same memory block if shared memory is activated. Input follows the DENSITY
 * enumeration from qmfunctions.h.
 */
void Density::alloc(int type) {
    switch (type) {
    case DENSITY::Total:
        if (hasTotal()) MSG_ERROR("Density not empty");
        if (isShared()) {
            this->dens_t = new FunctionTree<3>(*MRA, this->sh_mem);
        } else {
            this->dens_t = new FunctionTree<3>(*MRA);
        }
        break;
    case DENSITY::Spin:
        if (hasSpin()) MSG_ERROR("Density not empty");
        if (isShared()) {
            this->dens_s = new FunctionTree<3>(*MRA, this->sh_mem);
        } else {
            this->dens_s = new FunctionTree<3>(*MRA);
        }
        break;
    case DENSITY::Alpha:
        if (hasAlpha()) MSG_ERROR("Density not empty");
        if (isShared()) {
            this->dens_a = new FunctionTree<3>(*MRA, this->sh_mem);
        } else {
            this->dens_a = new FunctionTree<3>(*MRA);
        }
        break;
    case DENSITY::Beta:
        if (hasBeta()) MSG_ERROR("Density not empty");
        if (isShared()) {
            this->dens_b = new FunctionTree<3>(*MRA, this->sh_mem);
        } else {
            this->dens_b = new FunctionTree<3>(*MRA);
        }
        break;
    default:
        MSG_ERROR("Invalid density type");
        break;
    }
}

/** @brief clear density pointers
 *
 * @param type: density component
 *
 * This will set the FunctionTree pointer of the given density component to NULL
 * without deallocating. Input follows the DENSITY enumeration from qmfunctions.h.
 * Negative input clears all parts of the density.
 */
void Density::clear(int type) {
    if (type < 0 or type == DENSITY::Total) this->dens_t = 0;
    if (type < 0 or type == DENSITY::Spin)  this->dens_s = 0;
    if (type < 0 or type == DENSITY::Alpha) this->dens_a = 0;
    if (type < 0 or type == DENSITY::Beta)  this->dens_b = 0;
}

/** @brief free density pointers
 *
 * @param type: density component
 *
 * This will deallocate the FunctionTree pointer of the given density component.
 * Input follows the DENSITY enumeration from qmfunctions.h. Negative input frees
 * ALL components of the density.
 */
void Density::free(int type) {
    if (type < 0 or type == DENSITY::Total) { if (hasTotal()) delete this->dens_t; }
    if (type < 0 or type == DENSITY::Spin)  { if (hasSpin())  delete this->dens_s; }
    if (type < 0 or type == DENSITY::Alpha) { if (hasAlpha()) delete this->dens_a; }
    if (type < 0 or type == DENSITY::Beta)  { if (hasBeta())  delete this->dens_b; }
    clear(type);
}

/** @brief return number of MW nodes
 *
 * @param type: density component
 *
 * This will count the number of MW nodes of the given density component.
 * Input follows the DENSITY enumeration from qmfunctions.h. Negative input counts
 * ALL components of the density.
 */
int Density::getNNodes(int type) const {
    int nNodes = 0;
    if (type < 0 or type == DENSITY::Total) { if (hasTotal()) nNodes += total().getNNodes(); }
    if (type < 0 or type == DENSITY::Spin)  { if (hasSpin())  nNodes += spin().getNNodes(); }
    if (type < 0 or type == DENSITY::Alpha) { if (hasAlpha()) nNodes += alpha().getNNodes(); }
    if (type < 0 or type == DENSITY::Beta)  { if (hasBeta())  nNodes += beta().getNNodes(); }
    return nNodes;
}

/** @brief return meta data
 *
 * Tree sizes (nChunks) are flushed before return.
 */
DensityMeta& Density::getMetaData() {
    this->meta.nChunksTotal = 0;
    this->meta.nChunksSpin = 0;
    this->meta.nChunksAlpha = 0;
    this->meta.nChunksBeta = 0;
    if (hasTotal()) this->meta.nChunksTotal = total().getNChunksUsed();
    if (hasSpin())  this->meta.nChunksSpin  = spin().getNChunksUsed();
    if (hasAlpha()) this->meta.nChunksAlpha = alpha().getNChunksUsed();
    if (hasBeta())  this->meta.nChunksBeta  = beta().getNChunksUsed();
    return this->meta;
}

/** @brief set a particular density component function
 *
 * @param rho: density function
 * @param type: density component
 *
 * Input follows the DENSITY enumeration from qmfunctions.h. This will overwrite
 * any existing FunctionTree pointer.
 */
void Density::setDensity(FunctionTree<3> *rho, int type) {
         if (type == DENSITY::Total) this->dens_t = rho;
    else if (type == DENSITY::Spin)  this->dens_s = rho;
    else if (type == DENSITY::Alpha) this->dens_a = rho;
    else if (type == DENSITY::Beta)  this->dens_b = rho;
    else MSG_ERROR("Invalid density type");
}

} //namespace mrchem
