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
#include <string>

#include "MRCPP/Printer"

#include "AOBasis.h"
#include "Intgrl.h"
#include "OrbitalExp.h"
#include "chemistry/Nucleus.h"

using mrcpp::GaussExp;
using mrcpp::GaussFunc;
using mrcpp::Gaussian;

namespace mrchem {
namespace gto_utils {

OrbitalExp::OrbitalExp(Intgrl &intgrl)
        : cartesian(true) {
    readAOExpansion(intgrl);
}

OrbitalExp::~OrbitalExp() {
    for (auto &orbital : this->orbitals) {
        if (orbital != nullptr) { delete orbital; }
    }
}

GaussExp<3> OrbitalExp::getMO(int i, const DoubleMatrix &M) const {
    if (M.cols() != size()) MSG_ERROR("Size mismatch");
    GaussExp<3> mo_i;
    int n = 0;
    for (int j = 0; j < size(); j++) {
        GaussExp<3> ao_j = getAO(j);
        // ao_i.normalize();
        if (std::abs(M(i, j)) > mrcpp::MachineZero) {
            ao_j *= M(i, j);
            mo_i.append(ao_j);
            n++;
        }
    }
    if (n == 0) {
        MSG_WARN("No contributing orbital");
        GaussFunc<3> zeroFunc(0.0, 0.0);
        GaussExp<3> zeroExp;
        zeroExp.append(zeroFunc);
        mo_i.append(zeroExp);
    }
    // mo_i->normalize();
    return mo_i;
}

GaussExp<3> OrbitalExp::getDens(const DoubleMatrix &D) const {
    if (D.rows() != size()) MSG_ERROR("Size mismatch");
    if (D.cols() != size()) MSG_ERROR("Size mismatch");

    GaussExp<3> d_exp;
    for (int i = 0; i < size(); i++) {
        for (int j = 0; j < size(); j++) {
            GaussExp<3> ao_i = getAO(i);
            GaussExp<3> ao_j = getAO(j);
            GaussExp<3> d_ij = ao_i * ao_j;
            d_ij *= D(i, j);
            d_exp.append(d_ij);
        }
    }
    return d_exp;
}

void OrbitalExp::rotate(const DoubleMatrix &U) {
    std::vector<GaussExp<3> *> tmp;
    for (int i = 0; i < size(); i++) {
        auto *mo_i = new GaussExp<3>;
        int n = 0;
        for (int j = 0; j < size(); j++) {
            GaussExp<3> ao_j = getAO(j);
            // ao_j.normalize();
            if (std::abs(U(i, j)) > mrcpp::MachineZero) {
                ao_j *= U(i, j);
                mo_i->append(ao_j);
                n++;
            }
        }
        if (n == 0) {
            MSG_WARN("No contributing orbital");
            GaussFunc<3> zeroFunc(0.0, 0.0);
            GaussExp<3> ao_j;
            ao_j.append(zeroFunc);
            mo_i->append(ao_j);
        }
        // mo_i->normalize();
        tmp.push_back(mo_i);
    }
    for (int i = 0; i < size(); i++) {
        delete orbitals[i];
        orbitals[i] = tmp[i];
        tmp[i] = nullptr;
    }
}

void OrbitalExp::readAOExpansion(Intgrl &intgrl) {
    for (int i = 0; i < intgrl.getNNuclei(); i++) {
        Nucleus &nuc = intgrl.getNucleus(i);
        AOBasis &aoBasis = intgrl.getAOBasis(i);
        for (int j = 0; j < aoBasis.getNFunc(); j++) {
            GaussExp<3> *ao = new GaussExp<3>(aoBasis.getAO(j, nuc.getCoord()));
            this->orbitals.push_back(ao);
        }
    }
    transformToSpherical();
}

void OrbitalExp::transformToSpherical() {
    if (not this->cartesian) { return; }
    std::vector<GaussExp<3> *> tmp;
    int nOrbs = this->size();
    int n = 0;
    while (n < nOrbs) {
        int l = getAngularMomentum(n);
        if (l < 2) {
            GaussExp<3> *orb = this->orbitals[n];
            tmp.push_back(orb);
            this->orbitals[n] = nullptr;
            n++;
        } else if (l == 2) {
            int nprim = this->orbitals[n]->size();

            for (int i = 0; i < 6; i++) {
                if (this->orbitals[n + i]->size() != nprim) { MSG_ABORT("Contracted d orbials with different number of primitives"); }
            }

            auto *sph_xy = new GaussExp<3>;
            auto *sph_yz = new GaussExp<3>;
            auto *sph_z2 = new GaussExp<3>;
            auto *sph_xz = new GaussExp<3>;
            auto *sph_x2 = new GaussExp<3>;

            double sqrt3 = std::sqrt(3.0);

            for (int i = 0; i < nprim; i++) {
                Gaussian<3> &xx = this->orbitals[n + 0]->getFunc(i);
                Gaussian<3> &xy = this->orbitals[n + 1]->getFunc(i);
                Gaussian<3> &xz = this->orbitals[n + 2]->getFunc(i);
                Gaussian<3> &yy = this->orbitals[n + 3]->getFunc(i);
                Gaussian<3> &yz = this->orbitals[n + 4]->getFunc(i);
                Gaussian<3> &zz = this->orbitals[n + 5]->getFunc(i);

                sph_xy->append(xy);
                sph_yz->append(yz);

                sph_z2->append(xx);
                sph_z2->getFunc(sph_z2->size() - 1).setCoef(-0.5 * xx.getCoef());
                sph_z2->append(yy);
                sph_z2->getFunc(sph_z2->size() - 1).setCoef(-0.5 * yy.getCoef());
                sph_z2->append(zz);
                sph_z2->getFunc(sph_z2->size() - 1).setCoef(zz.getCoef());

                sph_xz->append(xz);

                sph_x2->append(xx);
                sph_x2->getFunc(sph_x2->size() - 1).setCoef(0.5 * sqrt3 * xx.getCoef());
                sph_x2->append(yy);
                sph_x2->getFunc(sph_x2->size() - 1).setCoef(-0.5 * sqrt3 * yy.getCoef());
            }

            tmp.push_back(sph_xy);
            tmp.push_back(sph_yz);
            tmp.push_back(sph_z2);
            tmp.push_back(sph_xz);
            tmp.push_back(sph_x2);

            n += 6;
        } else if (l == 3) {
            GaussExp<3> *sph[7];

            for (int i = 0; i < 7; i++) { sph[i] = new GaussExp<3>; }

            // order of cartesian f orbitals:
            // xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz

            double c1 = std::sqrt(2.5), c2 = std::sqrt(15.0), c3 = std::sqrt(1.5);

            // from page 211 of Molecular Electronic Structure Theory (Helgaker, et. al.)
            double coeffs[7][10] = {
                {0.0, 1.5 * c1, 0.0, 0.0, 0.0, 0.0, -0.5 * c1, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0, c2, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, -0.5 * c3, 0.0, 0.0, 0.0, 0.0, -0.5 * c3, 0.0, 2.0 * c3, 0.0},
                {0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0, -1.5, 0.0, 1.0},
                {-0.5 * c3, 0.0, 0.0, -0.5 * c3, 0.0, 2.0 * c3, 0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.5 * c2, 0.0, 0.0, 0.0, 0.0, -0.5 * c2, 0.0, 0.0},
                {0.5 * c1, 0.0, 0.0, -1.5 * c1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            };

            double normalization[10] = {15.0, 3.0, 3.0, 3.0, 1.0, 3.0, 15.0, 3.0, 3.0, 15.0};
            for (int i = 0; i < 10; i++) { normalization[i] = std::sqrt(normalization[i]); }

            int nprim = this->orbitals[n]->size();

            for (int i = 0; i < 10; i++) {
                if (this->orbitals[n + i]->size() != nprim) { MSG_ABORT("Contracted f orbials with different number of primitives"); }
            }

            for (int i = 0; i < nprim; i++) {
                for (int j = 0; j < 7; j++) {
                    for (int k = 0; k < 10; k++) {
                        if (coeffs[j][k] != 0.0) {
                            Gaussian<3> &func = this->orbitals[n + k]->getFunc(i);
                            sph[j]->append(func);
                            sph[j]->getFunc(sph[j]->size() - 1).setCoef(coeffs[j][k] * normalization[k] * func.getCoef());
                        }
                    }
                }
            }

            for (int i = 0; i < 7; i++) {
                sph[i]->normalize();
                std::cout << "Norm " << i << ": " << sph[i]->calcSquareNorm() << std::endl;
                tmp.push_back(sph[i]);
            }

            n += 10;
        } else {
            MSG_ABORT("Only s, p, d, and f orbitals are supported");
        }
    }
    for (int i = 0; i < nOrbs; i++) {
        if (this->orbitals[i] != nullptr) {
            delete this->orbitals[i];
            this->orbitals[i] = nullptr;
        }
    }
    this->orbitals.clear();
    for (auto &i : tmp) {
        this->orbitals.push_back(i);
        i = nullptr;
    }
    this->cartesian = false;
}

int OrbitalExp::getAngularMomentum(int n) const {
    int l = -1;
    GaussExp<3> &gExp = *this->orbitals[n];
    for (int i = 0; i < gExp.size(); i++) {
        const auto &pow = gExp.getPower(i);
        int iL = pow[0] + pow[1] + pow[2];
        if (l < 0) {
            l = iL;
        } else if (iL != l) {
            MSG_ABORT("Orbital is not pure angular momentum function");
        }
    }
    return l;
}

} // namespace gto_utils
} // namespace mrchem
