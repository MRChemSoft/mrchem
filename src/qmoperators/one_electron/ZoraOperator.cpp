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

#include "ZoraOperator.h"

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

ZoraOperator::ZoraOperator(QMPotential &vz, double c, double proj_prec, bool inverse) {
    Timer timer;
    double two_cc = 2.0 * c * c;

    auto k = std::make_shared<QMPotential>(1);
    qmfunction::deep_copy(*k, vz);

    if (k->hasImag()) MSG_ERROR("Inverse of complex function");
    if (k->hasReal()) {
        mrcpp::refine_grid(k->real(), 1);
        if (inverse) {
            k->real().map([two_cc](double val) { return (two_cc - val) / two_cc; });
        } else {
            k->real().map([two_cc](double val) { return two_cc / (two_cc - val); });
        }
        k->real().crop(proj_prec);
    }

    RankZeroOperator &kappa = (*this);
    kappa = k;
    if (inverse) {
        kappa.name() = "kappa_m1";
    } else {
        kappa.name() = "kappa";
    }
    auto plevel = Printer::getPrintLevel();
    print_utils::qmfunction(2, "ZORA operator (" + kappa.name() + ")", *k, timer);
}

} // namespace mrchem
