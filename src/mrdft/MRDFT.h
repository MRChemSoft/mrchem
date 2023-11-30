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

#include <memory>

#include <nlohmann/json.hpp>

#include "Functional.h"
#include "Grid.h"

namespace mrdft {

class MRDFT final {
public:
    MRDFT(std::unique_ptr<Grid> &g, std::unique_ptr<Functional> &f)
            : G(std::move(g))
            , F(std::move(f)) {}
    ~MRDFT() = default;

    auto &grid() { return *G; }
    auto &functional() { return *F; }

    mrcpp::FunctionTreeVector<3> evaluate(mrcpp::FunctionTreeVector<3> &inp);

private:
    std::unique_ptr<Grid> G{nullptr};
    std::unique_ptr<Functional> F{nullptr};
};

} // namespace mrdft
