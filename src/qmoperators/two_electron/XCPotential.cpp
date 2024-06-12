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

#include "XCPotential.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

using QMOperator_p = std::shared_ptr<mrchem::QMOperator>;

namespace mrchem {

/** @brief Prepare the operator for application
 *
 * @param[in] prec Apply precision
 *
 * Sequence of steps required to compute the XC potentials:
 *
 * 1) Compute density
 * 2) Setup xcfun input functions (gradients etc.)
 * 3) Evaluate xcfun
 * 4) Compute XC energy by integrating energy density
 * 5) Compute XC potential(s) from xcfun output functions
 *
 */
void XCPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    Timer timer;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(3, "Building XC operator");
    mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
    mrcpp::print::separator(3, '-');
    if (this->mrdft == nullptr) MSG_ERROR("XCFunctional not initialized");
    if (this->potentials.size() != 0) MSG_ERROR("Potential not properly cleared");

    auto &grid = this->mrdft->grid().get();
    mrcpp::FunctionTreeVector<3> xc_inp = setupDensities(prec, grid);
    this->xc_vec = this->mrdft->evaluate(xc_inp);

    // Fetch energy
    this->energy = this->mrdft->functional().XCenergy;

    // Fetch potential
    auto &v_local = mrcpp::get_func(this->xc_vec, 1);
    auto *v_global = new mrcpp::FunctionTree<3>(v_local.getMRA());
    mrcpp::copy_grid(*v_global, v_local);
    mrcpp::copy_func(*v_global, v_local);
    this->potentials.push_back(std::make_tuple(1.0, v_global));

    // Fetch potential
    if (this->mrdft->functional().isSpin()) {
        auto &v_local = mrcpp::get_func(this->xc_vec, 2);
        auto *v_global = new mrcpp::FunctionTree<3>(v_local.getMRA());
        mrcpp::copy_grid(*v_global, v_local);
        mrcpp::copy_func(*v_global, v_local);
        this->potentials.push_back(std::make_tuple(1.0, v_global));
    }


    if (plevel == 2) {
        int totNodes = 0;
        int totSize = 0;
        for (auto i = 0; i < this->potentials.size(); i++) {
            auto &f_i = mrcpp::get_func(this->potentials, i);
            totNodes += f_i.getNNodes();
            totSize += f_i.getSizeNodes();
        }
        auto t = timer.elapsed();
        mrcpp::print::tree(2, "XC operator", totNodes, totSize, t);
    }
    mrcpp::print::footer(3, timer, 2);
}

/** @brief Clears all data in the XCPotential object */
void XCPotential::clear() {
    this->energy = 0.0;
    for (auto &rho : this->densities) rho.free(NUMBER::Total);
    mrcpp::clear(this->potentials, true);
    mrcpp::clear(this->xc_vec, true);
    clearApplyPrec();
}

Density &XCPotential::getDensity(DensityType spin, int pert_idx) {
    int dens_idx = -1;
    if (spin == DensityType::Total) {
        if (pert_idx == 0) dens_idx = 0;
        if (pert_idx == 1) dens_idx = 3;
    } else if (spin == DensityType::Alpha) {
        if (pert_idx == 0) dens_idx = 1;
        if (pert_idx == 1) dens_idx = 4;
    } else if (spin == DensityType::Beta) {
        if (pert_idx == 0) dens_idx = 2;
        if (pert_idx == 1) dens_idx = 5;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }
    if (dens_idx < 0) MSG_ABORT("Invalid density index");
    if (dens_idx > densities.size()) MSG_ABORT("Invalid density index");
    return densities[dens_idx];
}

/** @brief Return FunctionTree for the XC spin potential
 *
 * @param[in] type Which spin potential to return (alpha, beta or total)
 */
FunctionTree<3> &XCPotential::getPotential(int spin) {
    int nPots = this->potentials.size();
    if (nPots < 1 or nPots > 2) MSG_ERROR("Invalid potential");

    bool spinFunctional = this->mrdft->functional().isSpin();
    int pot_idx = -1;
    if (spinFunctional and spin == SPIN::Alpha) {
        pot_idx = 0;
    } else if (spinFunctional and spin == SPIN::Beta) {
        pot_idx = 1;
    } else if (not spinFunctional) {
        pot_idx = 0;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }
    return mrcpp::get_func(this->potentials, pot_idx);
}

/** @brief XCPotentialD1 application
 *
 * @param[in] phi Orbital to which the potential is applied
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
 */
Orbital XCPotential::apply(Orbital phi) {
    QMPotential &V = *this;
    if (V.hasImag()) MSG_ERROR("Imaginary part of XC potential non-zero");

    FunctionTree<3> &pot = getPotential(phi.spin());
    V.setReal(&pot);
    Orbital Vphi = QMPotential::apply(phi);
    V.setReal(nullptr);
    return Vphi;
}

Orbital XCPotential::dagger(Orbital phi) {
    QMPotential &V = *this;
    if (V.hasImag()) MSG_ERROR("Imaginary part of XC potential non-zero");

    FunctionTree<3> &pot = getPotential(phi.spin());
    V.setReal(&pot);
    Orbital Vphi = QMPotential::dagger(phi);
    V.setReal(nullptr);
    return Vphi;
}

QMOperatorVector XCPotential::apply(QMOperator_p &O) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
