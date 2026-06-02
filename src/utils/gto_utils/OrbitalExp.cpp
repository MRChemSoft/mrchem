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
#include "AOContraction.h"
#include "Intgrl.h"
#include "OrbitalExp.h"
#include "chemistry/Nucleus.h"
#include "utils/math_utils.h"

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

GaussExp<3> OrbitalExp::getMO(int i, const DoubleMatrix &M, const double threshold) const {
    if (M.cols() != size()) MSG_ERROR("Size mismatch");
    GaussExp<3> mo_i;
    int n = 0;
    for (int j = 0; j < size(); j++) {
        GaussExp<3> ao_j = getAO(j);
        // ao_i.normalize();
        if (std::abs(M(i, j)) > threshold) {
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

/**
 * Computes the Ns coefficient from equation 6.4.49 on page 215 of the "Molecular Electronic Structure Theory" Book
 */
static double ns_coeff(int l, int m) {
    double lf = math_utils::factorial(l);
    double f1 = math_utils::factorial(l + m) / lf;
    double f2 = math_utils::factorial(l - m) / lf;
    double f3 = math_utils::pow_by_squaring(2.0, -(2 * std::abs(m) + (m == 0 ? 1 : 0) - 1));

    return std::sqrt(f1 * f2 * f3);
}

/**
 * Computes the C_tuv^lm coefficient from equation 6.4.48 on page 215 of the "Molecular Electronic Structure Theory" Book
 */
static double c_coeff(int l, int m, int t, int u, int v) {
    uint64_t c = math_utils::binomial(l, t) * math_utils::binomial(l - t, std::abs(m) + t) *
                 math_utils::binomial(t, u) *
                 math_utils::binomial(std::abs(m), 2 * v + (m < 0 ? 1 : 0));

    double cd = static_cast<double>(c);

    return ((t + v) % 2 == 0 ? cd : -cd) * math_utils::pow_by_squaring(4.0, -t);
}

/**
 * This computes directly the index a cartesian orbital given only the exponents of y and z.
 * (Note that the exponent of x is not needed)
 *
 * For example, the f orbitals (l = 3) come in the order:
 * cartesian_to_index(0, 0) == 0 (x^3)
 * cartesian_to_index(1, 0) == 1 (x^2 y)
 * cartesian_to_index(0, 1) == 2 (x^2 z)
 * cartesian_to_index(2, 0) == 3 (x y^2)
 * cartesian_to_index(1, 1) == 4 (x y z)
 * cartesian_to_index(0, 2) == 5 (x z^2)
 * cartesian_to_index(3, 0) == 6 (y^3)
 * cartesian_to_index(2, 1) == 7 (y^2 z)
 * cartesian_to_index(1, 2) == 8 (y z^2)
 * cartesian_to_index(0, 3) == 9 (z^3)
 *
 * @param ly exponent of y
 * @param lz exponent of z
 */
static int cartesian_to_index(int ly, int lz) {
    return (ly * (ly + 2 * lz + 1) + lz * (lz + 3)) / 2;
}

/**
 * Computes the coefficients for the solid spherical harmonics using
 * equation 6.4.47 on page 215 of the "Molecular Electronic Structure Theory" Book
 */
static CartToSphTransformation initializeSphCoeffs(int l) {
    CartToSphTransformation transformation;

    for (int m = -l; m <= l; m++) {
        std::vector<int> inds;
        std::vector<double> coeffs;

        int ml0 = m < 0 ? 1 : 0;
        double ns = ns_coeff(l, m);

        for (int t = 0; t <= (l - std::abs(m)) / 2; t++) {
            for (int u = 0; u <= t; u++) {
                for (int v = 0; v <= (std::abs(m) - ml0) / 2; v++) {
                    int lx = 2 * t + std::abs(m) - 2 * u - 2 * v - ml0;
                    int ly = 2 * u + 2 * v + ml0;
                    int lz = l - 2 * t - std::abs(m);

                    int ind = cartesian_to_index(ly, lz);

                    double coeff = c_coeff(l, m, t, u, v);
                    double norm = cartesianNormFac(lx, ly, lz);

                    inds.push_back(ind);
                    coeffs.push_back(ns * coeff * norm);
                }
            }
        }

        transformation.inds.push_back(inds);
        transformation.coeffs.push_back(coeffs);
    }

    return transformation;
}

CartToSphTransformation &OrbitalExp::getSphTransformation(int l) {
    while (sph_transformation_data.size() <= l) {
        int new_l = sph_transformation_data.size();
        CartToSphTransformation new_transformation = initializeSphCoeffs(new_l);
        sph_transformation_data.push_back(new_transformation);
    }

    return sph_transformation_data[l];
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
        } else {
            std::vector<GaussExp<3> *> sph;

            int ncart = ((l + 1) * (l + 2)) / 2;
            int nsph = 2 * l + 1;

            for (int i = 0; i < nsph; i++) { sph.push_back(new GaussExp<3>); }

            int nprim = this->orbitals[n]->size();

            CartToSphTransformation &trans_data = getSphTransformation(l);

            for (int i = 0; i < nprim; i++) {
                for (int j = 0; j < nsph; j++) {
                    std::vector<int> &inds = trans_data.inds[j];
                    std::vector<double> &coeffs = trans_data.coeffs[j];

                    for (int k = 0; k < inds.size(); k++) {
                        int ind = inds[k];
                        double coeff = coeffs[k];

                        Gaussian<3> &func = this->orbitals[n + ind]->getFunc(i);
                        sph[j]->append(func);
                        sph[j]->getFunc(sph[j]->size() - 1).setCoef(coeff * func.getCoef());
                    }
                }
            }

            for (int i = 0; i < nsph; i++) {
                tmp.push_back(sph[i]);
            }

            n += ncart;
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
    for (size_t i = 0; i < gExp.size(); i++) {
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
