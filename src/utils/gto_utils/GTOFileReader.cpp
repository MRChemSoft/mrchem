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
#include "qmfunctions/Orbital.h"

#ifdef MRCHEM_HAS_HDF5 
#include "utils/gto_utils/DiracAOReader.h"
#endif

#include <MRCPP/Gaussians>
#include <MRCPP/Parallel>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include <iomanip>
#include <sstream>

namespace mrchem {
namespace gto_utils {

GTOFileReader::GTOFileReader(const std::string &GTOType, const std::string &basisFile) {
    if (GTOType == "dalton")
        this->reader = std::make_unique<DaltonAOReader>(basisFile);
    else if (GTOType == "orca")
        this->reader = std::make_unique<DaltonAOReader>(basisFile);
    else if (GTOType == "dirac") {
#ifdef MRCHEM_HAS_HDF5
        this->reader = std::make_unique<DiracAOReader>(basisFile);
#else
        MSG_ABORT("MRChem was built without HDF5 support, which is required for gto_type = '" << GTOType
                  << "'. Reconfigure with --hdf5=default.");
#endif
    } else
        MSG_ABORT("GTOType " << GTOType << " not supported");
}

void GTOFileReader::projectMO(OrbitalVector &Phi, const std::string &mo_file, double prec, double screen, size_t nSpinors) {
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

    // Read MO file (transpose)
    mrcpp::Timer t1;
    std::vector<ComplexMatrix> MOCoefficients = this->reader->getMOCoefficients(mo_file, Phi.size(), nSpinors);
    std::vector<DoubleMatrix> Re(nSpinors), Im(nSpinors);
    for (size_t c = 0; c < nSpinors; c++) {
        Re[c] = MOCoefficients[c].real().transpose(); // nMOs x nAO
        Im[c] = MOCoefficients[c].imag().transpose();
    }
    OrbitalExp gto_exp(*this->reader);     // large (or only) component
    std::unique_ptr<OrbitalExp> small_exp; // small component, only for 4-component spinors
    if (MOCoefficients.size() == 4)
        small_exp = std::make_unique<OrbitalExp>(*this->reader, true);
    t1.stop();

    mrcpp::Timer t2;
    for (size_t i = 0; i < Phi.size(); i++) {
        mrcpp::Timer t_i;
        bool grid_was_built = false;
        if (mrcpp::mpi::my_func(Phi[i])) {
            if (nSpinors == 1) {
                grid_was_built = projectRealPart(Phi[i].real(), gto_exp, i, Re[0], prec, screen);
            } else {
                // DIRAC: nComp complex components
                Phi[i].free();       // remove any existing trees
                Phi[i].defcomplex(); // isreal = 0, iscomplex = 1
                for (size_t c = 0; c < nSpinors; c++) {
                    const OrbitalExp &exp = (c < 2) ? gto_exp : *small_exp; // components 0,1 large, 2,3 small
                    mrcpp::FunctionTree<3, double> reTree(*MRA);
                    mrcpp::FunctionTree<3, double> imTree(*MRA);
                    grid_was_built |= projectRealPart(reTree, exp, i, Re[c], prec, screen);
                    grid_was_built |= projectRealPart(imTree, exp, i, Im[c], prec, screen);
                    // re + i*im as one complex tree; the orbital takes ownership
                    Phi[i].setCplx(new mrcpp::FunctionTree<3, ComplexDouble>(reTree, imTree), c);
                }
            }
        }
        std::stringstream o_txt;
        o_txt << std::setw(w1 - 1) << i;
        o_txt << std::setw(w3) << print_utils::dbl_to_str(Phi[i].norm(), pprec, true);
        if (grid_was_built)
            o_txt << '*';
        print_utils::qmfunction(2, o_txt.str(), Phi[i], t_i);
    }
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    t2.stop();
    mrcpp::print::separator(2, '-');
    mrcpp::print::time(1, "Reading MO matrix", t1);
    mrcpp::print::time(1, "Projecting GTO MOs", t2);
}

bool GTOFileReader::projectRealPart(mrcpp::FunctionTree<3, double> &out, const OrbitalExp &exp, int i, const DoubleMatrix &M, double prec, double screen) {
    out.setZero();
    if (M.row(i).cwiseAbs().maxCoeff() < mrcpp::MachineZero)
        return false; // this part is exactly zero: keep the zero tree

    mrcpp::GaussExp<3> mo = exp.getMO(i, M);
    mo.calcScreening(screen);
    mrcpp::project(prec, out, mo);

    // Compare with the analytic norm of the Gaussian expansion (works for every component,
    // unlike the old check against 1, which only holds for a complete normalised orbital)
    const double refNorm = std::sqrt(mo.calcSquareNorm());
    if (std::abs(std::sqrt(out.getSquareNorm()) - refNorm) <= 1e2 * prec)
        return false;

    // Fallback, as before: refine the grid with a less screened expansion, then project again
    mrcpp::GaussExp<3> mo_grid = exp.getMO(i, M, prec * 1e-1);
    mrcpp::GaussExp<3> mo_screened = exp.getMO(i, M, prec * 1e-2);
    mrcpp::build_grid(out, mo_grid);
    mrcpp::project(prec, out, mo_screened);
    out.crop(prec * 1e-1, false);
    return true;
}

} // namespace gto_utils
} // namespace mrchem