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

#include "utils/gto_utils/GTOFileReader.h"

#include "utils/print_utils.h"
#include "utils/gto_utils/DaltonAOReader.h"
#include "utils/gto_utils/OrbitalExp.h"
#include "qmfunctions/Orbital.h"

#include <MRCPP/Gaussians>
#include <MRCPP/Parallel>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

namespace mrchem {
namespace gto_utils {

GTOFileReader::GTOFileReader(const std::string &GTOType, const std::string &basisFile) {
    if (GTOType == "Dalton")
        this->reader = std::make_unique<DaltonAOReader>(basisFile);
    else if (GTOType == "ORCA")
        this->reader = std::make_unique<DaltonAOReader>(basisFile);
    else    
        MSG_ABORT("GTOType " << GTOType << " not supported");
}

void GTOFileReader::projectMO(OrbitalVector &Phi, const std::string &mo_file, double prec, double screen) {
    auto pprec = mrcpp::Printer::getPrecision();
    auto w0 = mrcpp::Printer::getWidth() - 2;
    auto w1 = 5;
    auto w2 = w0 * 2 / 9;
    auto w3 = w0 - w1 - 3 * w2;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w3) << "Norm";
    o_head << std::setw(w2 + 1) << "Nodes";
    o_head << std::setw(w2) << "Size";
    o_head << std::setw(w2) << "Time";
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    // Setup AO basis
    OrbitalExp gto_exp(*this->reader);

    // Read MO file (transpose)
    mrcpp::Timer t1;
    DoubleMatrix MO = this->reader->getMOCoefficients(mo_file);
    if (static_cast<size_t>(MO.cols()) < Phi.size())
        MSG_ABORT("Size mismatch");
    t1.stop();

    mrcpp::Timer t2;
    for (size_t i = 0; i < Phi.size(); i++) {
        mrcpp::Timer t_i;
        bool grid_was_built = false;
        if (mrcpp::mpi::my_func(Phi[i])) {
            mrcpp::GaussExp<3> mo_i = gto_exp.getMO(i, MO.transpose());
            mo_i.calcScreening(screen);
            mrcpp::project(prec, Phi[i].real(), mo_i);
            if (std::abs(Phi[i].norm() - 1.0) > 1e2 * prec) {
                mrcpp::GaussExp<3> mo_i_grid = gto_exp.getMO(i, MO.transpose(), prec * 1e-1);
                mrcpp::GaussExp<3> mo_i_screened = gto_exp.getMO(i, MO.transpose(), prec * 1e-2);
                grid_was_built = true;
                mrcpp::build_grid(Phi[i].real(), mo_i_grid);
                mrcpp::project(prec, Phi[i].real(), mo_i_screened);
                Phi[i].real().crop(prec * 1e-1, false);
            }
        }
        std::stringstream o_txt;
        o_txt << std::setw(w1 - 1) << i;
        o_txt << std::setw(w3) << print_utils::dbl_to_str(Phi[i].norm(), pprec, true);
        if (grid_was_built)
            o_txt << '*';
        print_utils::qmfunction(2, o_txt.str(), Phi[i], t_i);
    }
    t2.stop();
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    mrcpp::print::separator(2, '-');
    mrcpp::print::time(1, "Reading MO matrix", t1);
    mrcpp::print::time(1, "Projecting GTO MOs", t2);
}

} // namespace gto_utils
} // namespace mrchem