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

#include <MRCPP/MWFunctions>

namespace mrdft {

class Grid final {
public:
    Grid(const mrcpp::MultiResolutionAnalysis<3> &mra)
            : tree(mra) {}
    ~Grid() {}

    auto &get() { return tree; }
    auto size() const { return tree.getNEndNodes(); }

    auto generate(int n) {
        mrcpp::FunctionTreeVector<3> out;
        for (int i = 0; i < n; i++) {
            auto *tmp = new mrcpp::FunctionTree<3>(tree.getMRA());
            mrcpp::copy_grid(*tmp, tree);
            out.push_back(std::make_tuple(1.0, tmp));
        }
        return out;
    }

    void unify(mrcpp::FunctionTreeVector<3> &inp) {
        // Extend current grid
        for (auto i = 0; i < inp.size(); i++) {
            auto &inp_i = mrcpp::get_func(inp, i);
            tree.appendTreeNoCoeff(inp_i);
        }
        // Unify input grids
        for (auto i = 0; i < inp.size(); i++) {
            auto &inp_i = mrcpp::get_func(inp, i);
            while (mrcpp::refine_grid(inp_i, tree)) {};
        }
    }

private:
    mrcpp::FunctionTree<3> tree;
};

} // namespace mrdft
