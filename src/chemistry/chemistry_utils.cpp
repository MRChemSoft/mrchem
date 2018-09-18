#include "chemistry_utils.h"
#include "Nucleus.h"
#include "utils/math_utils.h"

namespace mrchem {

/** @brief computes the repulsion self energy of a set of nuclei
 * 
 * @param[in] nucs the set of nuclei
 *
 */    
double compute_nuclear_repulsion(const Nuclei &nucs) {
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int i = 0; i < nNucs; i++) {
        const Nucleus &nuc_i = nucs[i];
        const double Z_i = nuc_i.getCharge();
        const double *R_i = nuc_i.getCoord();
        for (int j = i+1; j < nNucs; j++) {
            const Nucleus &nuc_j = nucs[j];
            const double Z_j = nuc_j.getCharge();
            const double *R_j = nuc_j.getCoord();
            double R_ij = math_utils::calc_distance(R_i, R_j);
            E_nuc += (Z_i*Z_j)/R_ij;
        }
    }
    return E_nuc;
}

} //namespace mrchem
