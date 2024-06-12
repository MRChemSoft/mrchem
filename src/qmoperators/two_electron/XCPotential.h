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

#include "qmoperators/QMPotential.h"
#include "qmfunctions/Density.h"

#include "mrdft/MRDFT.h"

/**
 * @class XCPotential
 * @brief Exchange-Correlation potential defined by a particular (spin) density
 *
 * The XC potential is computed by mapping of the density through a XC functional,
 * provided by the XCFun library. There are two ways of defining the density:
 *
 *  1) Use getDensity() prior to setup() and build the density as you like.
 *  2) Provide a default set of orbitals in the constructor that is used to
 *     compute the density on-the-fly in setup().
 *
 * If a set of orbitals has NOT been given in the constructor, the density
 * MUST be explicitly computed prior to setup(). The density will be computed
 * on-the-fly in setup() ONLY if it is not already available. After setup() the
 * operator will be fixed until clear(), which deletes both the density and the
 * potential.
 *
 * LDA and GGA functionals are supported as well as two different ways to compute
 * the XC potentials: either with explicit derivatives or gamma-type derivatives.
 *
 */

namespace mrchem {

class XCPotential : public QMPotential {
public:
    explicit XCPotential(std::unique_ptr<mrdft::MRDFT> &F, std::shared_ptr<OrbitalVector> Phi = nullptr, bool mpi_shared = false)
            : QMPotential(1, mpi_shared)
            , energy(0.0)
            , orbitals(Phi)
            , mrdft(std::move(F)) {}
    ~XCPotential() override = default;

    mrcpp::FunctionTreeVector<3> xc_vec;
    void setup(double prec) override;
    void clear() override;

    friend class XCOperator;

protected:
    double energy;                           ///< XC energy
    std::vector<Density> densities;          ///< XC densities (total or alpha/beta)
    mrcpp::FunctionTreeVector<3> potentials; ///< XC Potential functions collected in a vector
    std::shared_ptr<OrbitalVector> orbitals; ///< External set of orbitals used to build the density
    std::unique_ptr<mrdft::MRDFT> mrdft;     ///< External XC functional to be used

    double getEnergy() const { return this->energy; }
    Density &getDensity(DensityType spin, int pert_idx);
    mrcpp::FunctionTree<3> &getPotential(int spin);


    virtual mrcpp::FunctionTreeVector<3> setupDensities(double prec, mrcpp::FunctionTree<3> &grid) = 0;

    Orbital apply(Orbital phi) override;
    Orbital dagger(Orbital phi) override;
    QMOperatorVector apply(std::shared_ptr<QMOperator> &O) override;
};

} // namespace mrchem
