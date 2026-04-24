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

#pragma once

#include <nlohmann/json.hpp>

#include "SCFSolver.h"
#include "properties/SCFEnergy.h"
#include "external_solvers/ExternalSolver.h"
#include "external_solvers/ChemTensorSolver.h"


/** @class LagrangianSolver
 *
 * @brief Ground state SCF and multi-configurational optimization with Lagrangian optimization
 *
 * This ground state solver performs a constraint optimization of the orbitals
 * where the constraint is the orthonormality of the orbitals.
 * The solver utilizes an external driver to compute the 1-2RDMs.
 */

namespace mrchem {

class Molecule;
class FockBuilder;

class LagrangianSolver : public SCFSolver {
public:
    LagrangianSolver();
    virtual ~LagrangianSolver() override = default;

    
    OrbitalVector get_orbitals(){ return *this->orbitals; };

    nlohmann::json optimize(Molecule &mol, FockBuilder &F, ChemTensorSolver &S);

protected:
    //std::string chkFile;  ///< Name of checkpoint file
    int nIter{};
    float scf_tol{};
    std::vector<double> energy{};
    std::shared_ptr<OrbitalVector> orbitals{};
    

    void set_orbitals(OrbitalVector Phi);

    void orbital_update(ChemTensorSolver &S);
    std::shared_ptr<ComplexMatrix> calculate_lagrange_multipliers(ChemTensorSolver &S);

};

} // namespace mrchem




/*  Open questions:
    - Store orbitals separately or override the ones in molecule?
    - save solver as a private element?
    
*/