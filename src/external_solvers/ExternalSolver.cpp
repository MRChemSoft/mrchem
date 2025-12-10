/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/MWOperators"

#include "ExternalSolver.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/GenericTwoOrbitalsOperator.h"

namespace mrchem {

class FockBuilder;
// class PoissonOperator;

/** @brief Calculates and stores the one- and two-electron integrals for given orbitals
 *
 * @param Phi: Vector of orbitals
 *
 * Calculates the one- and two-electron integrals for the orbitals in Phi, and stores them
 * in the class members one_body_integrals and two_body_integrals.
 *
 */

// TODO: more efficient without defining new operators?
void ExternalSolver::set_integrals(OrbitalVector &Phi, FockBuilder &F) {
    F.setup(this->prec);
    // operators
    KineticOperator K(F.momentum());
    NuclearOperator V = *(F.getNuclearOperator());
    GenericTwoOrbitalsOperator g = *(F.getGenericTwoOrbitalsOperator());

    g.setup(std::make_shared<OrbitalVector>(Phi), this->prec);
    // set the one- and two-body integrals
    ExternalSolver::set_one_body_integrals(Phi, K, V);
    ExternalSolver::set_two_body_integrals(Phi, g);
}

// Private

// TODO: change 'NuclearOperator' to 'RankZeroOperator'
void ExternalSolver::set_one_body_integrals(OrbitalVector &Phi, KineticOperator &K, NuclearOperator &V) {
    this->one_body_integrals = std::make_shared<ComplexMatrix>(K(Phi, Phi) + V(Phi, Phi));
}

void ExternalSolver::set_two_body_integrals(OrbitalVector &Phi, GenericTwoOrbitalsOperator &g) {
    int n_orb = Phi.size();
    this->two_body_integrals = std::make_shared<ComplexTensorR4>(n_orb, n_orb, n_orb, n_orb);
    this->two_body_integrals->setZero();
    std::cout << "two body integral size: " << this->two_body_integrals->dimensions() << std::endl;
    // TODO: use 8-fold symmetry
    for (int j = 0; j < n_orb; j++) {
        for (int l = 0; l < n_orb; l++) {
            g.set_pair(j, l);
            for (int k = 0; k < n_orb; k++) {
                // calculate |g_jl|Phi_k>
                // Orbital tmp_k = Phi[k].paramCopy();
                // mrcpp::cplxfunc::multiply(tmp_i, Phi[i].dagger(), Vjl, this->prec, true, true);
                Orbital tmp_k = g.apply(Phi[k]);
                for (int i = 0; i < n_orb; i++) {
                    // calculate (ij|kl) = <Phi_i|V_jl|Phi_k>
                    // BUG: add factor 4pi??
                    (*this->two_body_integrals)(i, j, k, l) = mrcpp::dot(Phi[i], tmp_k);
                }
            }
        }
    }
}

} // namespace mrchem