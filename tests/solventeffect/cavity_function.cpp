/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "catch.hpp"

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Plotter"
#include <array>
#include <cmath>
#include <vector>

#include "mrchem.h"
#include "parallel.h"

#include "chemistry/Cavity.h"

using namespace mrchem;

namespace cavity_function {

TEST_CASE("Cavityfunction", "[cavity_function]") {
    // const double prec = 1.0e-3;
    const double prec = 1.0e-3;
    const double thrs = 1.0e-8;

    // test volume of single sphere
    std::vector<mrcpp::Coord<3>> coords = {{0.0, 0.0, 0.0}};
    std::vector<double> R = {1.0};
    double slope = 0.2;

    Cavity sphere(coords, R, slope);
    mrcpp::FunctionTree<3> cav_tree(*MRA);
    mrcpp::project<3>(prec, cav_tree, sphere);

    auto sphere_volume = cav_tree.integrate();
    REQUIRE(sphere_volume == Approx(4.4401176759).epsilon(thrs));

    // test gradCavity
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    auto dcav = sphere.getGradVector();
    mrcpp::FunctionTree<3> dx_cav_analytical(*MRA);
    mrcpp::project<3>(prec, dx_cav_analytical, dcav[0]);

    mrcpp::FunctionTreeVector<3> d_cav = mrcpp::gradient(*D_p, cav_tree);

    mrcpp::FunctionTree<3> diff_func(*MRA);
    mrcpp::add(-1.0, diff_func, 1.0, mrcpp::get_func(d_cav, 0), -1.0, dx_cav_analytical);
    std::cout << diff_func.integrate() << std::endl;
    int aPts = 100;                     // Number of points in a
    int bPts = 100;                     // Number of points in b
    mrcpp::Coord<3> o{-2.0, -2.0, 0.0}; // Origin vector
    mrcpp::Coord<3> a{4.0, 0.0, 0.0};   // Boundary vector A
    mrcpp::Coord<3> b{0.0, 4.0, 0.0};
    mrcpp::Plotter<3> plot(o);                                           // Plotter of 3D functions
    plot.setRange(a, b);                                                 // Set plot range
    plot.surfPlot({aPts, bPts}, mrcpp::get_func(d_cav, 0), "dx_cav");    // Write to file f_tree.surf
    plot.surfPlot({aPts, bPts}, dx_cav_analytical, "dx_cav_analytical"); // Write to file f_tree.surf
    plot.surfPlot({aPts, bPts}, diff_func, "diff_func");                 // Write to file f_tree.surf

    mrcpp::Coord<3> lo{-2.0, 0.0, 0.0}; // Origin vector
    mrcpp::Plotter<3> line(lo);
    line.setRange(a);
    line.linePlot({aPts}, cav_tree, "cav_line");
    line.linePlot({aPts}, mrcpp::get_func(d_cav, 0), "dx_cav_line");
    line.linePlot({aPts}, dx_cav_analytical, "dx_cav_analytical_line");
    line.linePlot({aPts}, diff_func, "diff_func_line");

    // test volume of two interlocking spheres
    coords.push_back({0.0, 0.0, 1.0});
    R.push_back(1.0);
    Cavity two_spheres(coords, R, slope);
    mrcpp::FunctionTree<3> two_cav_tree(*MRA);
    mrcpp::project<3>(prec, two_cav_tree, two_spheres);

    auto two_sphere_volume = two_cav_tree.integrate();
    REQUIRE(two_sphere_volume == Approx(7.50966)); //.epsilon(thrs))
}
} // namespace cavity_function
