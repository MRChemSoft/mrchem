/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "periodic_utils.h"

#include "MRCPP/Printer"

#include "chemistry/Nucleus.h"
#include "utils/math_utils.h"

namespace mrchem {

std::array<double, 3> periodic::on_boundary(double period, mrcpp::Coord<3> coord) {
    std::array<double, 3> counter{0.0};
    for (auto i = 0; i < 3; i++) {
        if (abs(coord[i] - period * 0.5) < 1.0e-4) counter[i] = -1.0;
        if (abs(coord[i] + period * 0.5) < 1.0e-4) counter[i] = 1.0;
    }
    return counter;
}

Nuclei periodic::fix_boundary_charge(Nuclei nucs, double period) {

    for (auto &nuc : nucs) {
        auto boundary_vector = periodic::on_boundary(period, nuc.getCoord());
        auto boundary_count = static_cast<int>(math_utils::calc_distance(boundary_vector, {}));
        auto Z = nuc.getCharge();
        if (boundary_count == 1) nuc.setCharge(Z / 2.0);
        if (boundary_count == 2) nuc.setCharge(Z / 4.0);
        if (boundary_count == 3) nuc.setCharge(Z / 8.0);
    }

    return nucs;
}

Nuclei periodic::periodify_nuclei(Nuclei nucs, double period, double shift) {

    Nuclei new_nucs;
    auto top = period / 2.0;
    for (auto &nuc : nucs) {
        auto element = nuc.getElement().getSymbol();
        auto coord = nuc.getCoord();
        auto ob = periodic::on_boundary(period, coord);
        std::vector<mrcpp::Coord<3>> new_coords;
        for (auto x = 0; x < 2; x++) {
            for (auto y = 0; y < 2; y++) {
                for (auto z = 0; z < 2; z++) {
                    coord = nuc.getCoord();
                    coord[0] += ob[0] * x * period;
                    coord[1] += ob[1] * y * period;
                    coord[2] += ob[2] * z * period;
                    for (auto i = 0; i < 3; i++) {
                        if (ob[i] != 0.0) coord[i] *= ob[i] * ob[i] * shift;
                    }
                    new_coords.push_back(coord);
                }
            }
        }
        sort(new_coords.begin(), new_coords.end());
        new_coords.erase(unique(new_coords.begin(), new_coords.end()), new_coords.end());
        for (auto &c : new_coords) new_nucs.push_back(element, c);
    }
    return new_nucs;
}

double periodic::calc_rc(Nuclei nucs, double period) {
    std::vector<mrcpp::Coord<3>> new_coords;
    for (auto &nuc : nucs) {
        for (auto x = -1; x < 2; x++) {
            for (auto y = -1; y < 2; y++) {
                for (auto z = -1; z < 2; z++) {
                    auto coord = nuc.getCoord();
                    coord[0] += x * period;
                    coord[1] += y * period;
                    coord[2] += z * period;
                    new_coords.push_back(coord);
                }
            }
        }
    }
    std::vector<double> dists;
    for (auto c_i : new_coords) {
        for (auto b_i : new_coords) {
            auto dist = math_utils::calc_distance(c_i, b_i);
            if (dist > 0.00001) dists.push_back(dist);
        }
    }
    return *min_element(dists.begin(), dists.end()) * 0.499999;
}

} // namespace mrchem
