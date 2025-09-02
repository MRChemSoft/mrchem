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

#include <iostream>
#include <vector>

#include "MRCPP/MWFunctions"
#include <MRCPP/Printer>

#include "Element.h"
#include "PeriodicTable.h"
#include "pseudopotential/pseudopotential.h"

namespace mrchem {

class Nucleus final {
public:
    Nucleus(const Element &elm, const mrcpp::Coord<3> &r, double rms = -1.0)
            : charge(elm.getZ())
            , radius(rms)
            , coord(r)
            , element(&elm) {}

    Nucleus(const Element &elm, double z, const mrcpp::Coord<3> &r,double rms = -1.0)
            : charge(z)
            , radius(rms)
            , coord(r)
            , element(&elm) {}

    Nucleus(const Element &elm, const mrcpp::Coord<3> &r, std::shared_ptr<PseudopotentialData> pp_data, double rms = -1.0) {
        charge = pp_data->getZeff();
        radius = rms;
        coord = r;
        element = &elm;
        pp_data_ptr = pp_data;
        has_pp_data = true;
    }
    
    
    Nucleus(const Nucleus &nuc)
            : charge(nuc.charge)
            , radius(nuc.radius)
            , coord(nuc.coord)
            , element(nuc.element)
            , pp_data_ptr(nuc.pp_data_ptr)
            , has_pp_data(nuc.has_pp_data) {}
    
    Nucleus &operator=(const Nucleus &nuc) {
        if (this != &nuc) {
            this->charge = nuc.charge;
            this->radius = nuc.radius;
            this->coord = nuc.coord;
            this->element = nuc.element;
            this->pp_data_ptr = nuc.pp_data_ptr;
            this->has_pp_data = nuc.has_pp_data;
        }
        return *this;
    }

    /**
     * Set the charge of the nucleus. Deprecated.
     */
    void setCharge(double z) {
        // throw std::runtime_error("Nucleus::setCharge() is deprecated");
        MSG_ABORT("Nucleus::setCharge() is deprecated why should the charge be changed?");
        this->charge = z;
        }

    void setRMSRadius(double r) {
        this->radius = r; }

    void setCoord(const mrcpp::Coord<3> &r) { this->coord = r; }

    /**
     * Get the charge of the nucleus. Returns the effective charge if a pseudopotential is used.
     */
    double getCharge() const { return this->charge; }

    /**
     * get the atomic number of the nucleus. When pseudopotentials are used the output is different from getCharge()
     */
    double getAtomicNumber() const { return this->element->getZ(); }

    /**
     * Get the pseudopotential data of the nucleus.
     */
    std::shared_ptr<PseudopotentialData> getPseudopotentialData() const {
        if (!has_pp_data) {
            MSG_ABORT("Nucleus has no pseudopotential data");
        }
        return this->pp_data_ptr;
    }

    /**
     * Check if the nucleus has pseudopotential data.
     */
    bool hasPseudopotential() const { return this->has_pp_data; }


    double getRMSRadius() const { return this->radius; }
    const mrcpp::Coord<3> &getCoord() const { return this->coord; }
    const Element &getElement() const { return *this->element; }

    friend std::ostream &operator<<(std::ostream &o, const Nucleus &nuc) {
        o << std::endl << *nuc.element << "   Z = ";
        o << nuc.charge << std::endl << "   ";
        for (int d = 0; d < 3; d++) o << nuc.coord[d] << " ";
        o << std::endl;
        return o;
    }

private:
    double charge;
    double radius;
    mrcpp::Coord<3> coord;
    const Element *element;
    std::shared_ptr<PseudopotentialData> pp_data_ptr;
    bool has_pp_data = false;


};

class Nuclei : public std::vector<Nucleus> {
public:
    void push_back(const Nucleus &nuc) { std::vector<Nucleus>::push_back(nuc); }
    void push_back(const std::string &atom, const mrcpp::Coord<3> &xyz, double rms = -1.0) {
        PeriodicTable pt;
        Nucleus nuc(pt.getElement(atom.c_str()), xyz, rms);
        std::vector<Nucleus>::push_back(nuc);
    }
    void push_back(const std::string &atom, const mrcpp::Coord<3> &xyz, std::shared_ptr<PseudopotentialData> pp_data, double rms = -1.0) {
        PeriodicTable pt;
        Nucleus nuc(pt.getElement(atom.c_str()), xyz, pp_data, rms);
        std::vector<Nucleus>::push_back(nuc);
    }
};

} // namespace mrchem
