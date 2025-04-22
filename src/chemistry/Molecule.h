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
#include "properties/GeometricDerivative.h"
#include "properties/Magnetizability.h"
#include "properties/NMRShielding.h"
#include "properties/OrbitalEnergies.h"
#include "properties/Polarizability.h"
#include "properties/QuadrupoleMoment.h"
#include "properties/SCFEnergy.h"
#include "qmfunctions/Orbital.h"
#include "properties/HirshfeldCharges.h"

/** @class Molecule
 *
 * @brief Collection of properties related to a molecule.
 *
 * A molecule is basically a collection of nuclei and a collection of properties.
 * Mainly used for bookkeeping and printing of properties.
 *
 */

namespace mrchem {
class Cavity;
template <typename P> using PropertyMap = std::map<std::string, P>;

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

    mrcpp::Coord<3> calcCenterOfMass() const;
    mrcpp::Coord<3> calcCenterOfCharge() const;

    /**
     * @brief get all nuclei
     */
    auto &getNuclei() { return this->nuclei; }

    /**
     * @brief get_pseudo_potential_nuclei
     */
    Nuclei getPseudoPotentialNuclei();


    /**
     * @brief return true if the molecule has a pseudopotential, false otherwise
     */
    bool hasPseudopotential() const;

    /**
     * @brief return true if the molecule has a NLCC pseudopotential, false otherwise
     */
    bool hasNLCCPseudopotential() const;

    /**
     * @brief return trun if the molecule has a pseudopotential with projectors, false otherwise
     */
    bool hasProjectorPseudopotential() const;

    /**
     * @brief get all electron nuclei
     */
    Nuclei getAllElectronNuclei();


    auto &getOrbitals() { return *this->orbitals_0; }
    auto &getOrbitalsX() { return *this->orbitals_x; }
    auto &getOrbitalsY() { return *this->orbitals_y; }
    auto &getCavity() { return *this->cavity; }
    auto &getFockMatrix() { return this->fock_matrix; }

    const auto &getNuclei() const { return this->nuclei; }
    const auto &getOrbitals() const { return *this->orbitals_0; }
    const auto &getOrbitalsX() const { return *this->orbitals_x; }
    const auto &getOrbitalsY() const { return *this->orbitals_y; }
    const auto &getFockMatrix() const { return this->fock_matrix; }

    auto getOrbitals_p() const { return this->orbitals_0; }
    auto getOrbitalsX_p() const { return this->orbitals_x; }
    auto getOrbitalsY_p() const { return this->orbitals_y; }
    auto getCavity_p() const { return this->cavity; }

    nlohmann::json json() const;
    void printGeometry() const;
    void printEnergies(const std::string &txt) const;
    void printProperties() const;

    void initPerturbedOrbitals(bool dynamic);
    void initCavity(const std::vector<mrcpp::Coord<3>> &coords, const std::vector<double> &R, const std::vector<double> &alphas, const std::vector<double> &betas, const std::vector<double> &sigmas);

    SCFEnergy &getSCFEnergy() { return this->energy; }
    OrbitalEnergies &getOrbitalEnergies() { return this->epsilon; }
    DipoleMoment &getDipoleMoment(const std::string &id) { return this->dipole.at(id); }
    QuadrupoleMoment &getQuadrupoleMoment(const std::string &id) { return this->quadrupole.at(id); }
    Polarizability &getPolarizability(const std::string &id) { return this->polarizability.at(id); }
    Magnetizability &getMagnetizability(const std::string &id) { return this->magnetizability.at(id); }
    NMRShielding &getNMRShielding(const std::string &id) { return this->nmr_shielding.at(id); }
    GeometricDerivative &getGeometricDerivative(const std::string &id) { return this->geometric_derivative.at(id); }
    HirshfeldCharges &getHirshfeldCharges(const std::string &id) { return this->hirshfeld_charges.at(id); }

    PropertyMap<DipoleMoment> &getDipoleMoments() { return this->dipole; }
    PropertyMap<QuadrupoleMoment> &getQuadrupoleMoments() { return this->quadrupole; }
    PropertyMap<Polarizability> &getPolarizabilities() { return this->polarizability; }
    PropertyMap<Magnetizability> &getMagnetizabilities() { return this->magnetizability; }
    PropertyMap<NMRShielding> &getNMRShieldings() { return this->nmr_shielding; }
    PropertyMap<GeometricDerivative> &getGeometricDerivatives() { return this->geometric_derivative; }
    PropertyMap<HirshfeldCharges> &getHirshfeldCharges() { return this->hirshfeld_charges; }

protected:
    int charge{0};
    int multiplicity{1};
    Nuclei nuclei{};
    ComplexMatrix fock_matrix{};

    std::shared_ptr<Cavity> cavity{nullptr};
    std::shared_ptr<OrbitalVector> orbitals_0{std::make_shared<OrbitalVector>()};
    std::shared_ptr<OrbitalVector> orbitals_x{nullptr};
    std::shared_ptr<OrbitalVector> orbitals_y{nullptr};

    // Properties
    SCFEnergy energy{};
    OrbitalEnergies epsilon{};
    PropertyMap<DipoleMoment> dipole{};
    PropertyMap<QuadrupoleMoment> quadrupole{};
    PropertyMap<Polarizability> polarizability{};
    PropertyMap<Magnetizability> magnetizability{};
    PropertyMap<NMRShielding> nmr_shielding{};
    PropertyMap<GeometricDerivative> geometric_derivative{};
    PropertyMap<HirshfeldCharges> hirshfeld_charges{};

    void readCoordinateFile(const std::string &file);
    void readCoordinateString(const std::vector<std::string> &coord_str);
};

} // namespace mrchem
