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

#include <fstream>

#include "Molecule.h"
#include "Nucleus.h"

#include "environment/Cavity.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::Coord;
using mrcpp::Printer;

namespace mrchem {

/** @brief Constructor
 *
 * @param c: total charge
 * @param m: spin multiplicity
 */
Molecule::Molecule(int c, int m)
        : charge(c)
        , multiplicity(m) {}

/** @brief Constructor
 *
 * @param coord_file: xyz file with nuclear coordinates
 * @param c: total charge
 * @param m: spin multiplicity
 *
 * Nuclei are copied, all properties are uninitialized at this point.
 */
Molecule::Molecule(const std::string &coord_file, int c, int m)
        : charge(c)
        , multiplicity(m) {
    readCoordinateFile(coord_file);
}

/** @brief Constructor
 *
 * @param coord_file: list of stings with nuclear coordinates
 * @param c: total charge
 * @param m: spin multiplicity
 *
 * Nuclei are copied, all properties are uninitialized at this point.
 */
Molecule::Molecule(const std::vector<std::string> &coord_str, int c, int m)
        : charge(c)
        , multiplicity(m) {
    readCoordinateString(coord_str);
}

void Molecule::initPerturbedOrbitals(bool dynamic) {
    if (dynamic) {
        this->orbitals_x = std::make_shared<OrbitalVector>();
        this->orbitals_y = std::make_shared<OrbitalVector>();
    } else {
        this->orbitals_x = std::make_shared<OrbitalVector>();
        this->orbitals_y = this->orbitals_x;
    }
}

/** @brief Return number of electrons */
int Molecule::getNElectrons() const {
    auto totZ = 0;
    for (auto i = 0; i < getNNuclei(); i++) totZ += getNuclei()[i].getElement().getZ();
    return totZ - this->charge;
}

/** @brief Compute nuclear center of mass */
Coord<3> Molecule::calcCenterOfMass() const {
    Coord<3> COM;
    COM.fill(0.0);

    auto M = 0.0;
    for (auto i = 0; i < getNNuclei(); i++) {
        const auto &nuc = getNuclei()[i];
        const auto &r_i = nuc.getCoord();
        const auto m_i = nuc.getElement().getMass();
        for (auto d = 0; d < 3; d++) COM[d] += r_i[d] * m_i;
        M += m_i;
    }
    for (auto d = 0; d < 3; d++) COM[d] *= 1.0 / M;
    return COM;
}

/** @brief Compute nuclear center of charge */
Coord<3> Molecule::calcCenterOfCharge() const {
    Coord<3> COC;
    COC.fill(0.0);

    auto Z = 0.0;
    for (auto i = 0; i < getNNuclei(); i++) {
        const auto &nuc = getNuclei()[i];
        const auto &r_i = nuc.getCoord();
        const auto z_i = nuc.getElement().getZ();
        for (auto d = 0; d < 3; d++) COC[d] += r_i[d] * z_i;
        Z += z_i;
    }
    for (auto d = 0; d < 3; d++) COC[d] *= 1.0 / Z;
    return COC;
}

/** @brief Read nuclear coordinates from xyz file
 *
 * First entry in file is number of atoms:
 *
 * nAtoms
 * symbol   x_coord     y_coord     z_coord
 * symbol   x_coord     y_coord     z_coord
 * symbol   x_coord     y_coord     z_coord
 *
 */
void Molecule::readCoordinateFile(const std::string &coord_file) {
    std::ifstream ifs(coord_file.c_str());
    if (not ifs) MSG_ABORT("Failed to open coordinate file: " << coord_file);

    int nNuclei;
    Coord<3> coord;
    std::string sym;
    ifs >> nNuclei;
    for (auto i = 0; i < nNuclei; i++) {
        ifs >> sym;
        ifs >> coord[0];
        ifs >> coord[1];
        ifs >> coord[2];
        this->nuclei.push_back(sym.c_str(), coord);
    }
    ifs.close();
}

/** @brief Read nuclear coordinates from vector of strings
 *
 * Each entry in the vector of strings contains one atom:
 *
 *      "symbol   x_coord     y_coord     z_coord"
 *
 */
void Molecule::readCoordinateString(const std::vector<std::string> &coord_str) {
    Coord<3> coord;
    std::string sym;
    for (const auto &i : coord_str) {
        std::stringstream ss;
        ss.str(i);
        ss >> sym;
        ss >> coord[0];
        ss >> coord[1];
        ss >> coord[2];
        this->nuclei.push_back(sym.c_str(), coord);
    }
}

/** @brief Pretty output of molecular geometry */
void Molecule::printGeometry() const {
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 8;
    auto w3 = 2 * w0 / 9;
    auto w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::string(w4 - 1, ' ') << ':';
    o_head << std::setw(w3) << "x";
    o_head << std::setw(w3) << "y";
    o_head << std::setw(w3) << "z";

    mrcpp::print::header(0, "Molecule");
    print_utils::scalar(0, "Charge", getCharge(), "", 0);
    print_utils::scalar(0, "Multiplicity", getMultiplicity(), "", 0);
    mrcpp::print::separator(0, '-');
    println(0, o_head.str());
    mrcpp::print::separator(0, '-');
    for (auto i = 0; i < getNNuclei(); i++) {
        const auto &nuc = getNuclei()[i];
        std::stringstream o_sym;
        o_sym << std::setw(w1 - 1) << i;
        o_sym << std::setw(w2) << nuc.getElement().getSymbol();
        print_utils::coord(0, o_sym.str(), nuc.getCoord());
    }

    mrcpp::print::separator(0, '-');
    print_utils::coord(0, "Center of mass", calcCenterOfMass());
    mrcpp::print::separator(0, '=', 2);
}

void Molecule::printEnergies(const std::string &txt) const {
    energy.print(txt);
    epsilon.print(txt);
}

/** @brief Pretty output of molecular properties
 *
 * Only properties that have been initialized will be printed.
 */
void Molecule::printProperties() const {
    // For each std::map entry (first = id string, second = property)
    // this will print the property with its id appearing in the header
    for (const auto &dip : dipole) dip.second.print(dip.first);
    for (const auto &qua : quadrupole) qua.second.print(qua.first);
    for (const auto &pol : polarizability) pol.second.print(pol.first);
    for (const auto &mag : magnetizability) mag.second.print(mag.first);
    for (const auto &nmr : nmr_shielding) nmr.second.print(nmr.first);
    for (const auto &geo : geometric_derivative) geo.second.print(geo.first);
    for (const auto &hir : hirshfeld_charges) hir.second.print(hir.first);
}

nlohmann::json Molecule::json() const {
    nlohmann::json json_out;
    json_out["geometry"] = {};
    for (auto i = 0; i < getNNuclei(); i++) {
        const auto &nuc = getNuclei()[i];
        nlohmann::json json_atom;
        json_atom["symbol"] = nuc.getElement().getSymbol();
        json_atom["xyz"] = nuc.getCoord();
        json_out["geometry"].push_back(json_atom);
    }
    json_out["charge"] = getCharge();
    json_out["multiplicity"] = getMultiplicity();
    json_out["center_of_mass"] = calcCenterOfMass();

    json_out["scf_energy"] = energy.json();
    json_out["orbital_energies"] = epsilon.json();
    if (not dipole.empty()) json_out["dipole_moment"] = {};
    if (not quadrupole.empty()) json_out["quadrupole_moment"] = {};
    if (not polarizability.empty()) json_out["polarizability"] = {};
    if (not magnetizability.empty()) json_out["magnetizability"] = {};
    if (not nmr_shielding.empty()) json_out["nmr_shielding"] = {};
    if (not geometric_derivative.empty()) json_out["geometric_derivative"] = {};
    for (const auto &dip : dipole) json_out["dipole_moment"][dip.first] = dip.second.json();
    for (const auto &qua : quadrupole) json_out["quadrupole_moment"][qua.first] = qua.second.json();
    for (const auto &pol : polarizability) json_out["polarizability"][pol.first] = pol.second.json();
    for (const auto &mag : magnetizability) json_out["magnetizability"][mag.first] = mag.second.json();
    for (const auto &nmr : nmr_shielding) json_out["nmr_shielding"][nmr.first] = nmr.second.json();
    for (const auto &geo : geometric_derivative) json_out["geometric_derivative"][geo.first] = geo.second.json();
    for (const auto &hir : hirshfeld_charges) json_out["hirshfeld_charges"][hir.first] = hir.second.json();

    return json_out;
}

void Molecule::initCavity(const std::vector<mrcpp::Coord<3>> &coords,
                          const std::vector<double> &R,
                          const std::vector<double> &alphas,
                          const std::vector<double> &betas,
                          const std::vector<double> &sigmas) {
    if (cavity) MSG_ABORT("Cavity already initialized");
    this->cavity = std::make_shared<Cavity>(coords, R, alphas, betas, sigmas);
}

} // namespace mrchem
