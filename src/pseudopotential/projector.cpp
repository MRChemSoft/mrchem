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

#include "pseudopotential/projector.h"
#include <math.h>
#include <fstream>
#include <MRCPP/Printer>
#include <MRCPP/functions/GaussFunc.h>

#include <string>

ProjectorFunction::ProjectorFunction(mrcpp::Coord<3> pos, double rl, int i, int l, int m, double prec) {
    this->pos = pos;
    this->rl = rl;
    this->i = i;
    int ii = i + 1;
    this->l = l;
    this->m = m;
    this->prec = prec;
    // select the spherical harmonic function based on the angular momentum and magnetic quantum number
    switch_sperics(l, m);
    double prefactor = std::sqrt(2.0) / (std::pow(rl, l + (4.0 * ii - 1) / 2.0) * std::sqrt(std::tgamma( l + (4.0 * ii - 1.0) / 2.0 )) );

    auto project_analytic = [this, prefactor, ii](const std::array<double, 3> &r) -> double {
        std::array<double, 3> rprime = {r[0] - this->pos[0], r[1] - this->pos[1], r[2] - this->pos[2]};
        double normr = std::sqrt( rprime[0] * rprime[0] + rprime[1] * rprime[1] + rprime[2] * rprime[2]);
        return prefactor * std::pow(normr, 2 * (ii - 1)) * std::exp(- 0.5 * normr * normr / (this->rl * this-> rl) ) * this->s(rprime, normr);
    };
    // auto op = (*this);
    // mrcpp::ComplexFunction f;
    projector_ptr = std::make_shared<mrcpp::CompFunction<3>>();

    double sigma = 0.6;
    double beta = 0.5 / (sigma * sigma);
    mrcpp::GaussFunc<3> gauss(beta, 1.0, this->pos);

    mrcpp::project(*projector_ptr, gauss, prec);
    mrcpp::project(*projector_ptr,  static_cast<std::function<double(const mrcpp::Coord<3>&)>>(project_analytic), prec);

    double nrm = projector_ptr->norm();


    if (std::abs(nrm - 1.0) > 10 * prec) {
        MSG_ABORT("Projection of projector failed: norm=" << nrm << ", nodes=" << projector_ptr->getNNodes() << ", l=" << l << ", m=" << m << ", rl=" << rl << ", i=" << i);
    }

}

void ProjectorFunction::switch_sperics(int l, int m){
    this->s = get_spherical_harmonics(l, m);
}
