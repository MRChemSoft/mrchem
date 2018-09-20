/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

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
