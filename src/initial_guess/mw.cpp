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

#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "mw.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace mw {

bool project_mo(OrbitalVector &Phi, double prec, const std::string &mo_file,  int &file_ncomp);

} // namespace mw
} // namespace initial_guess

bool initial_guess::mw::setup(OrbitalVector &Phi, double prec, const std::string &file_p, const std::string &file_a, const std::string &file_b, int n_components) {
    if (Phi.size() == 0) return false;
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation   ", "Compute initial orbitals");
    print_utils::text(0, "Method        ", "Project MW molecular orbitals");
    print_utils::text(0, "Precision     ", print_utils::dbl_to_str(prec, 5, true));
    if (! orbital::size_doubly(Phi)) {
        print_utils::text(0, "Restricted    ", "False");
        print_utils::text(0, "MO alpha file ", file_a);
        print_utils::text(0, "MO beta file  ", file_b);
    } else {
        print_utils::text(0, "Restricted    ", "True");
        print_utils::text(0, "MO file ", file_p);
    }
    mrcpp::print::separator(0, '~', 2);

    // Separate alpha/beta from paired orbitals
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);

    // Project paired, alpha and beta separately
    auto success = true;
    int file_ncomp_p = 1, file_ncomp_a = 1, file_ncomp_b = 1; //max number of components for each type of orbitals in the stored files
    success &= initial_guess::mw::project_mo(Phi, prec, file_p, file_ncomp_p);
    success &= initial_guess::mw::project_mo(Phi_a, prec, file_a, file_ncomp_a);
    success &= initial_guess::mw::project_mo(Phi_b, prec, file_b, file_ncomp_b);
    // For 2C: beta orbitals are Kramers partners of alpha; move their
    // content from comp[0] to comp[1] to match Kramers symmetry
    // If the loaded files for beta orbitals are scalar (1 component)
    // swap them over to the second component to follow symmetry
    // and emulate the time-reversal operator -iσ_y K0
    if (n_components > 1 && file_ncomp_b == 1) {
        for (auto &phi : Phi_b) {
            if (not mrcpp::mpi::my_func(phi)) continue;
            //swap components
            std::swap(phi.CompD[0], phi.CompD[1]);
            std::swap(phi.CompC[0], phi.CompC[1]);
            //swap metadata
            std::swap(phi.func_ptr->data.Nchunks[0], phi.func_ptr->data.Nchunks[1]);
            //multiplying prefactors with -i σ_y
            std::swap(phi.func_ptr->data.c1[0], phi.func_ptr->data.c1[1]);
            phi.func_ptr->data.c1[1] *= -1.0;
        }
    }

    // Collect orbitals into one vector
    Phi = orbital::adjoin(Phi, Phi_a);
    Phi = orbital::adjoin(Phi, Phi_b);

    return success;
}

bool initial_guess::mw::project_mo(OrbitalVector &Phi, double prec, const std::string &mo_file,  int &file_ncomp) {
    if (Phi.size() == 0) return true;

    Timer t_tot;
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 2;
    auto w1 = 5;
    auto w2 = w0 * 2 / 9;
    auto w3 = w0 - w1 - 3 * w2;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w3) << "Norm";
    o_head << std::setw(w2 + 1) << "Nodes";
    o_head << std::setw(w2) << "Size";
    o_head << std::setw(w2) << "Time";

    mrcpp::print::header(1, "MW Initial Guess");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    bool success = true;
    bool iscomplex = false;
    file_ncomp = 1; //tracking the number of components from the loaded files (default will be scalar)
    for (int i = 0; i < Phi.size(); i++) {
        Timer t_i;
        if (mrcpp::mpi::my_func(Phi[i])) {
            std::stringstream orbname;
            orbname << mo_file << "_idx_" << i;

            Orbital phi_i;
            orbital::loadOrbital(orbname.str(), phi_i);
            file_ncomp = std::max(file_ncomp, phi_i.Ncomp()); //Keeping the highest number (in practice all files should have the same number of components)
            if (phi_i.getSquareNorm() < 0.0) {
                MSG_ERROR("Guess orbital not found: " << orbname.str());
                success &= false;
            }
            //determine the number of components to load
            //Useful if trying to load scalar orbitals into 2C spinors
            int n_load = std::min(Phi[i].Ncomp(), phi_i.Ncomp()); 
            if (phi_i.isreal()) {
                Phi[i].defreal();
                for (int comp = 0; comp<Phi[i].Ncomp(); comp++) {
                    if (comp < n_load){
                        // Refine to get accurate function values
                        mrcpp::refine_grid(phi_i.real(comp), 1);
                        mrcpp::project(prec, Phi[i].real(comp), phi_i.real(comp));
                    } else Phi[i].alloc_comp(comp, true);
                }
            }
            if (phi_i.iscomplex()) {
                Phi[i].defcomplex();
                for (int comp = 0; comp<Phi[i].Ncomp(); comp++) {
                    if (comp < n_load){
                        // Refine to get accurate function values
                        mrcpp::refine_grid(phi_i.complex(comp), 1);
                        iscomplex = true;
                        mrcpp::project(prec, Phi[i].complex(comp), phi_i.complex(comp));
                    } else Phi[i].alloc_comp(comp, true);
                }
            }
            std::stringstream o_txt;
            o_txt << std::setw(w1 - 1) << i;
            o_txt << std::setw(w3) << print_utils::dbl_to_str(Phi[i].norm(), pprec, true);
            print_utils::qmfunction(1, o_txt.str(), Phi[i], t_i);
        }
    }
    // if any orbital of "master" is complex, all are set complex (also the one which are not "mine")
#ifdef MRCHEM_HAS_MPI
    MPI_Bcast(&iscomplex, 1, MPI_C_BOOL, 0, mrcpp::mpi::comm_wrk);
#endif
    // NOTE: I don't know much about MPI, so this fix was done using Claude.
    // Each rank only loads orbitals where my_func(Phi[i]) is true. In an MPI run, rank 0 might own zero beta orbitals — its file_ncomp stays 1 even if the files are 2C. But the Kramers swap in setup() runs on all ranks (no my_func guard, same as core/sad). So without a reduction, ranks that loaded nothing would incorrectly trigger the swap while others don't.
    // MPI_Allreduce with MPI_MAX gives all ranks the true maximum component count found across all loaded files. Bcast from rank 0 would be wrong here (rank 0 might have loaded nothing).
    // WARNING: Claude mentions that "The existing MPI_Bcast(&iscomplex, ...) above code has exactly this flaw — if rank 0 owns no orbitals, it broadcasts false even if rank 1 found complex orbitals. "
// #ifdef MRCHEM_HAS_MPI
//     int global_file_ncomp = file_ncomp;
//     MPI_Allreduce(&file_ncomp, &global_file_ncomp, 1, MPI_INT, MPI_MAX, mrcpp::mpi::comm_wrk);
//     file_ncomp = global_file_ncomp;
// #endif
    if (iscomplex) {
        for (int i = 0; i < Phi.size(); i++) Phi[i].defcomplex();
    }
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    mrcpp::print::footer(1, t_tot, 2);
    return success;
}

} // namespace mrchem
