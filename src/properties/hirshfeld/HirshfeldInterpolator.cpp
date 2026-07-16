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

#include "HirshfeldInterpolator.h"
#include "qmfunctions/density_utils.h"

double HirshfeldRadInterpolater::getNorm() const {
    double norm = 0.0;
    double xmax = 20;
    int n = 1000;
    double dx = xmax / n;
    for (int i = 0; i < n; i++) {
        double x = i * dx;
        double y = std::exp(evalf(x));
        norm += y * x * x;
    }
    norm *= dx * 4 * M_PI;
    return norm;
}

// Constructor
HirshfeldRadInterpolater::HirshfeldRadInterpolater(const std::string element, std::string data_dir, bool writeToFile) {
    Eigen::VectorXd rGrid;
    Eigen::VectorXd rhoGrid;

    std::string filename = data_dir + '/' + element + ".density";

    mrchem::density::readAtomicDensity(filename, rGrid, rhoGrid);

    rhoGrid = rhoGrid.array().log();

    lnRho = std::make_shared<interpolation_utils::PolyInterpolator>(rGrid, rhoGrid);
    if (writeToFile) { writeInterpolatedDensity(element + ".interpolated"); }
}

// Function to evaluate the interpolated function
double HirshfeldRadInterpolater::evalf(const double &r) const {
    return lnRho->evalfLeftNoRightLinear(r);
}

void HirshfeldRadInterpolater::writeInterpolatedDensity(const std::string path) {
    std::ofstream file;
    file.open(path);
    for (double r = 0.0; r < 50.0; r += 0.01) {
        double rho = this->evalf(r);
        file << r << " " << rho << std::endl;
    }
    file.close();
}