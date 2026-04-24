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

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "LagrangianSolver.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/two_electron/FockBuilder.h"
#include <iostream>

using mrcpp::Printer;
//using mrcpp::Timer;
using nlohmann::json;

namespace mrchem {

// TODO: read all the settings in input
LagrangianSolver::LagrangianSolver(){
    this->nIter = 1;
    this->scf_tol = 1e-3;
}

void LagrangianSolver::set_orbitals(OrbitalVector Phi_n){
    this->orbitals = std::make_shared<OrbitalVector>(Phi_n);
}


/** @brief Run Lagrangian orbital optimization
 *
 * @param mol: Molecule to optimize
 * @param F: Fock Builder
 * @param S: External solver to compute 1-2RDMs
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm performs
 * a Lagrangian optimization of the orbitals, where the constraint is the orthonormal
 * condition of the orbitals. The driver is used to compute the 1-2RDMs.
 */
json LagrangianSolver::optimize(Molecule &mol, FockBuilder &F, ChemTensorSolver &S) {
    //printParameters("Optimize ground state orbitals");

    //Timer t_tot;
    json json_out;

    //Nuclei nucs = mol.getNuclei();

    // initial orbitals
    set_orbitals(mol.getOrbitals());
    double prec = 1e-3;

    bool converged = false;
    
    for(int i=0; i<this->nIter && !converged; i++){
        // set integrals, run DMRG and calculate RDMS
        S.set_integrals(*this->orbitals);
        S.optimize();
        this->energy.push_back(S.get_energy());

        // update orbitals
        orbital_update();

        // check if converged
        if(i>0){
            if(abs(this->energy[i-1]-this->energy[i])<this->scf_tol)
                converged=true;
        }
    }

    // end
    F.clear();
    json_out["converged"] = converged;
    return json_out;

}

void LagrangianSolver::orbital_update(){
    return;
}




} // namespace mrchem