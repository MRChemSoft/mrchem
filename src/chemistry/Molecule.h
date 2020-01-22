/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

/**
 *
 * \date Jul 13, 2010
 * \author Stig Rune Jensen \n
 *		   CTCC, University of Tromsø
 *
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <MRCPP/Printer>
#include <nlohmann/json.hpp>

#include "Nucleus.h"
#include "properties/DipoleMoment.h"
#include "properties/Magnetizability.h"
#include "properties/NMRShielding.h"
#include "properties/Polarizability.h"
#include "properties/QuadrupoleMoment.h"
#include "properties/SCFEnergy.h"
#include "qmfunctions/Orbital.h"

/** @class Molecule
 *
 * @brief Collection of properties related to a molecule.
 *
 * A molecule is basically a collection of nuclei and a collection of properties.
 * Mainly used for bookkeeping and printing of properties.
 *
 */

namespace mrchem {

class Molecule final {
public:
    explicit Molecule(int c = 0, int m = 1);
    explicit Molecule(const std::string &coord_file, int c = 0, int m = 1);
    explicit Molecule(const std::vector<std::string> &coord_str, int c = 0, int m = 1);
    Molecule(const Molecule &mol) = delete;
    Molecule &operator=(const Molecule &mol) = delete;

    int getCharge() const { return this->charge; }
    int getMultiplicity() const { return this->multiplicity; }
    int getNNuclei() const { return this->nuclei.size(); }
    int getNElectrons() const;

    void setCharge(int c) { this->charge = c; }
    void setMultiplicity(int m) { this->multiplicity = m; }
    void setGaugeOrigin(const mrcpp::Coord<3> &o) { this->origin = o; }

    mrcpp::Coord<3> getGaugeOrigin() const { return this->origin; }
    mrcpp::Coord<3> calcCenterOfMass() const;
    mrcpp::Coord<3> calcCenterOfCharge() const;

    auto &getNuclei() { return this->nuclei; }
    auto &getOrbitals() { return *this->orbitals_0; }
    auto &getOrbitalsX() { return *this->orbitals_x; }
    auto &getOrbitalsY() { return *this->orbitals_y; }
    auto &getFockMatrix() { return this->fock_matrix; }

    const auto &getNuclei() const { return this->nuclei; }
    const auto &getOrbitals() const { return *this->orbitals_0; }
    const auto &getOrbitalsX() const { return *this->orbitals_x; }
    const auto &getOrbitalsY() const { return *this->orbitals_y; }
    const auto &getFockMatrix() const { return this->fock_matrix; }

    auto getOrbitals_p() const { return this->orbitals_0; }
    auto getOrbitalsX_p() const { return this->orbitals_x; }
    auto getOrbitalsY_p() const { return this->orbitals_y; }

    nlohmann::json json() const;
    void printGeometry() const;
    void printProperties() const;

    void initPerturbedOrbitals(bool dynamic);

    SCFEnergy &getSCFEnergy();
    DipoleMoment &getDipoleMoment();
    QuadrupoleMoment &getQuadrupoleMoment();
    Magnetizability &getMagnetizability();
    NMRShielding &getNMRShielding(int k);
    Polarizability &getPolarizability(double omega);

protected:
    int charge{0};
    int multiplicity{1};
    Nuclei nuclei{};
    mrcpp::Coord<3> origin{};
    ComplexMatrix fock_matrix{};

    std::shared_ptr<OrbitalVector> orbitals_0{std::make_shared<OrbitalVector>()};
    std::shared_ptr<OrbitalVector> orbitals_x{nullptr};
    std::shared_ptr<OrbitalVector> orbitals_y{nullptr};

    // Properties
    std::unique_ptr<SCFEnergy> energy{nullptr};
    std::unique_ptr<DipoleMoment> dipole{nullptr};
    std::unique_ptr<QuadrupoleMoment> quadrupole{nullptr};
    std::unique_ptr<Magnetizability> magnetizability{nullptr};
    std::vector<std::unique_ptr<Polarizability>> polarizability{};
    std::vector<std::unique_ptr<NMRShielding>> nmr{};

    void initNuclearProperties(int nNucs);
    void readCoordinateFile(const std::string &file);
    void readCoordinateString(const std::vector<std::string> &coord_str);
};

} // namespace mrchem
