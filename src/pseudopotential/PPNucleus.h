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

#include "analyticfunctions/NuclearFunction.h"
#include "chemistry/Nucleus.h"
#include "MRCPP/Printer"

#include "utils/math_utils.h"
// #include "pseudopotential/HGH.hpp"
#include "pseudopotential/pseudopotential.h"

namespace mrchem {

/**
 * @class PPNucleus
 * @brief Represents a pseudopotential nucleus. (local part of the pseudopotential)
 * 
 * The PPNucleus class is a derived class of the NuclearFunction class. It represents a pseudopotential nucleus
 * and provides methods for evaluating the nuclear function and calculating parameters.
 */
class PPNucleus : public NuclearFunction {
public:

    /**
     * @brief Constructs a PPNucleus object with the given pseudopotential data.
     * @param pps The pseudopotential data. (One for each nucleus)
     */
    PPNucleus(const Nuclei &nucs) : NuclearFunction(){
        for(int i = 0; i < nucs.size(); i++) {
            if (! nucs[i].hasPseudopotential()) {
                MSG_ABORT("Nucleus has no pseudopotential data in constroctur of PPNucleus");
            }
            this->pps.push_back(*nucs[i].getPseudopotentialData());
        }
    }

    /**
     * @brief Evaluates the localized potential at the given coordinate.
     * @param r The coordinate at which to evaluate the nuclear function.
     * @return The value of the localized function at the given coordinate.
     */
    double evalf(const mrcpp::Coord<3> &r) const override {
        double result = 0.0;
        double temp_exp;
        for (int i = 0; i < this->nuclei.size(); i++) {
            const auto &R = this->nuclei[i].getCoord();
            auto R1 = math_utils::calc_distance(R, r);

            temp_exp = std::exp(-this->pps[i].alpha_pp * this->pps[i].alpha_pp * R1 * R1);

            result -= this->pps[i].getZeff() / R1 * std::erf(this->pps[i].alpha_pp * R1);

            double temp = 1.0;
            double temp_square = 2.0 * R1 * R1 * this->pps[i].alpha_pp * this->pps[i].alpha_pp;
            
            for (int iloc = 0; iloc < pps[i].nloc; iloc++){
                result += this->pps[i].c[iloc] * temp_exp * temp;
                temp *= temp_square;
            }
        }
        return result;
    }

    // I had to implement this functions, thez return not meaningful values.
    std::string getParamName1() const { return "Precision"; }
    std::string getParamName2() const { return "Smoothing"; }
    double calcParam1(double prec, const Nucleus &nuc) const { return prec; }
    double calcParam2(double prec, const Nucleus &nuc) const {
        double tmp = 0.0;
        return prec;
    }

protected:
    std::vector<PseudopotentialData> pps;

};

} // namespace mrchem
