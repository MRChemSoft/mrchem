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

#include "MRCPP/utils/parallel.h"
#include <MRCPP/Printer>
#include <XCFun/xcfun.h>
#include <fstream>
#include <iostream>
#include "MRCPP/utils/parallel.h"

#include "mrchem.h"
#include "mrenv.h"
#include "utils/print_utils.h"
#include "utils/json_utils.h"   // robust, tolerant JSON helpers
#include "version.h"

// #ifdef MRCHEM_ENABLE_LIBXC
// #  include "LibXCBackend.h"
// #  include <xc.h>
// #endif

using json    = nlohmann::json;
using Printer = mrcpp::Printer;

namespace mrchem {

namespace mrenv {
void init_printer(const json &json_print);
void init_mra(const json &json_mra);
void init_mpi(const json &json_mpi);
void print_header();
} // namespace mrenv

// -----------------------------------------------------------------------------
// Read the program JSON and return the normalized "input" subtree
// -----------------------------------------------------------------------------
json mrenv::fetch_json(int argc, char **argv) {
    const char *infile = nullptr;
    if (argc == 1) {
        infile = "STDIN";
    } else if (argc == 2) {
        infile = argv[1];
    } else {
        MSG_ERROR("Invalid number of arguments! Pass either no argument (read from stdin) or a single JSON file.");
    }

    json root;
    if (std::string(infile) == "STDIN") {
        std::cin >> root;
    } else {
        std::ifstream ifs(infile, std::ios_base::in);
        if (!ifs) MSG_ABORT("Could not open input file: " << infile);
        ifs >> root;
        ifs.close();
    }

    if (!root.contains("input") || !root["input"].is_object()) {
        MSG_ABORT(
            "Malformed program JSON: missing object \"input\".\n"
            "Note: mrchem.x expects the generated program JSON (with \"input\" and \"output\" sections),\n"
            "not the user .inp/.json you write by hand. Use `mrchem --dryrun <name>` or\n"
            "`mrchem --json --stdout <name> > <name>.json` to generate the program JSON."
        );
    }

    // Work on a copy of the "input" subtree and normalize boolean-ish fields (0/1, on/off, etc.)
    json j = root["input"];
    mrchem::json_utils::sanitize_booleans(j);
    return j;
}

// -----------------------------------------------------------------------------
// Initialization (printer, MRA, MPI) with tolerant JSON reads
// -----------------------------------------------------------------------------
void mrenv::initialize(const json &json_inp) {
    auto json_print = json_inp.find("printer");
    auto json_mra   = json_inp.find("mra");
    auto json_mpi   = json_inp.find("mpi");

    if (json_mra == json_inp.end() || !json_mra->is_object()) {
        MSG_ABORT("Missing MRA input! The \"input\" JSON must contain an \"mra\" object.");
    } else {
        mrenv::init_mra(*json_mra);
    }
    if (json_mpi   != json_inp.end() && json_mpi->is_object())   mrenv::init_mpi(*json_mpi);
    if (json_print != json_inp.end() && json_print->is_object()) mrenv::init_printer(*json_print);

    mrenv::print_header();
}

void mrenv::init_printer(const json &json_print) {
    // Safe, typed reads with defaults
    int  print_level = json_utils::value_loose<int>(json_print,  "print_level", 1);
    int  print_prec  = json_utils::value_loose<int>(json_print,  "print_prec",  10);
    int  print_width = json_utils::value_loose<int>(json_print,  "print_width", 120);
    bool print_mpi   = json_utils::value_loose<bool>(json_print, "print_mpi",   false);

    std::string fname = "mrchem.out";
    if (json_print.contains("file_name") && json_print["file_name"].is_string()) {
        fname = json_print["file_name"].get<std::string>();
    }

    if (print_mpi) {
        Printer::init(print_level, mrcpp::mpi::world_rank, mrcpp::mpi::world_size, fname.c_str());
    } else {
        Printer::init(print_level, mrcpp::mpi::world_rank, mrcpp::mpi::world_size);
    }
    Printer::setPrecision(print_prec);
    Printer::setWidth(print_width);
}

void mrenv::init_mra(const json &json_mra) {
    // Required fields: enforce exact types here (schema-generated JSON should be correct)
    int min_scale = json_mra.at("min_scale").get<int>();
    int max_scale = json_mra.at("max_scale").get<int>();
    auto corner   = json_mra.at("corner");
    auto boxes    = json_mra.at("boxes");
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes);

    int         order = json_mra.at("basis_order").get<int>();
    std::string btype = json_mra.at("basis_type").get<std::string>();

    auto max_depth = max_scale - min_scale;
    if (min_scale < mrcpp::MinScale) MSG_ABORT("Root scale too large");
    if (max_scale > mrcpp::MaxScale) MSG_ABORT("Max scale too large");
    if (max_depth > mrcpp::MaxDepth) MSG_ABORT("Max depth too large");

    if (btype == "interpolating") {
        mrcpp::InterpolatingBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else if (btype == "legendre") {
        mrcpp::LegendreBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else {
        MSG_ABORT("Invalid basis type!");
    }
    mrcpp::SetdefaultMRA(MRA);
}

void mrenv::init_mpi(const json &json_mpi) {
    mrcpp::mpi::numerically_exact = json_mpi["numerically_exact"];
    mrcpp::mpi::shared_memory_size = json_mpi["shared_memory_size"];
    mrcpp::mpi::bank_size = json_mpi["bank_size"];
    mrcpp::mpi::bank_per_node = json_mpi["bank_per_node"];
    mrcpp::mpi::omp_threads = json_mpi["omp_threads"];
    mrcpp::mpi::use_omp_num_threads = json_mpi["use_omp_num_threads"];
    mrcpp::mpi::initialize(); // NB: must be after bank_size and init_mra but before init_printer and print_header
}

// -----------------------------------------------------------------------------
// Pretty header
// -----------------------------------------------------------------------------
void mrenv::print_header() {
    auto pwidth      = Printer::getWidth();
    auto txt_width   = 50;
    auto pre_spaces  = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str  = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');
    std::stringstream o_ver, o_branch, o_hash, o_author, o_date;
    o_ver    << "VERSION            " << program_version();
    o_branch << "Git branch         " << git_branch();
    o_hash   << "Git commit hash    " << git_commit_hash();
    o_author << "Git commit author  " << git_commit_author();
    o_date   << "Git commit date    " << git_commit_date();

    int ver_len    = static_cast<int>(o_ver.str().size());
    int branch_len = static_cast<int>(o_branch.str().size());
    int hash_len   = static_cast<int>(o_hash.str().size());
    int auth_len   = static_cast<int>(o_author.str().size());
    int date_len   = static_cast<int>(o_date.str().size());

    o_ver    << std::string(std::max(0, txt_width - ver_len),    ' ');
    o_branch << std::string(std::max(0, txt_width - branch_len), ' ');
    o_hash   << std::string(std::max(0, txt_width - hash_len),   ' ');
    o_author << std::string(std::max(0, txt_width - auth_len),   ' ');
    o_date   << std::string(std::max(0, txt_width - date_len),   ' ');

    std::stringstream o_bank;
    if (mrcpp::mpi::bank_size > 0) {
        o_bank << "(" << mrcpp::mpi::tot_bank_size << " bank)";
    } else {
        o_bank << "(no bank)";
    }

    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '*');
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << " __  __ ____   ____ _                             " << post_str);
    println(0, pre_str << "|  \\/  |  _ \\ / ___| |__   ___ _ __ ___           " << post_str);
    println(0, pre_str << "| |\\/| | |_) | |   | '_ \\ / _ \\ '_ ` _ \\          " << post_str);
    println(0, pre_str << "| |  | |  _ <| |___| | | |  __/ | | | | |         " << post_str);
    println(0, pre_str << "|_|  |_|_| \\_\\\\____|_| |_|\\___|_| |_| |_|         " << post_str);
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << o_ver.str() << post_str);
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << o_branch.str() << post_str);
    println(0, pre_str << o_hash.str() << post_str);
    println(0, pre_str << o_author.str() << post_str);
    println(0, pre_str << o_date.str() << post_str);
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << "Contact: luca.frediani@uit.no                     " << post_str);
    println(0, pre_str << "                                                  " << post_str);
    println(0, pre_str << "Radovan Bast            Magnar Bjorgve            " << post_str);
    println(0, pre_str << "Roberto Di Remigio      Antoine Durdek            " << post_str);
    println(0, pre_str << "Luca Frediani           Gabriel Gerez             " << post_str);
    println(0, pre_str << "Stig Rune Jensen        Jonas Juselius            " << post_str);
    println(0, pre_str << "Rune Monstad            Peter Wind                " << post_str);
    println(0, pre_str << "                                                  " << post_str);
    mrcpp::print::separator(0, '*', 1);
    mrcpp::print::separator(0, '-', 1);
    print_utils::scalar(0, "MPI processes  ", mrcpp::mpi::world_size, o_bank.str(), 0, false);
    print_utils::scalar(0, "OpenMP threads ", mrcpp::omp::n_threads, "", 0, false);
    print_utils::scalar(0, "Total cores    ",
                        (mrcpp::mpi::world_size - mrcpp::mpi::tot_bank_size) * mrcpp::omp::n_threads
                        + mrcpp::mpi::tot_bank_size, "", 0, false);
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '-', 1);
    //     #ifdef MRCHEM_ENABLE_LIBXC
    // if (mrcpp::mpi::grand_master()) {
    //     println(0, "LibXC version " << XC_VERSION << " (" << XC_RELEASE_DATE << ")");
    //     println(0, "LibXC compiled with " << xc_compiled_functionals() << " functionals");
    // }
    // else printout(0, xcfun_splash());
    // #endif
    printout(0, xcfun_splash());
    mrcpp::print::environment(0);
    MRA->print();
}

void mrenv::finalize(double wt) {
    if (MRA != nullptr) delete MRA;
    MRA = nullptr;

    auto pwidth      = Printer::getWidth();
    auto txt_width   = 45;
    auto pre_spaces  = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str  = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');

    auto hr  = static_cast<int>(wt / 3600.0);
    auto min = static_cast<int>(std::fmod(wt, 3600.0) / 60.0);
    auto sec = static_cast<int>(std::fmod(wt, 60.0));

    std::stringstream o_time;
    o_time << "Wall time : " << std::setw(2) << hr << "h" << std::setw(3) << min << "m" << std::setw(3) << sec << "s";

    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '*');
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << "                Exiting MRChem               " << post_str);
    println(0, pre_str << "                                             " << post_str);
    println(0, pre_str << "           " << o_time.str() << "           " << post_str);
    println(0, pre_str << "                                             " << post_str);
    mrcpp::print::separator(0, '*');
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, ' ');
}

void mrenv::dump_json(const json &json_inp, const json &json_out) {
    json json_tot;
    json_tot["input"]  = json_inp;
    json_tot["output"] = json_out;

    // Try to derive a base name from printer.file_name; fallback to "mrchem"
    std::string base = "mrchem";
    try {
        if (json_inp.contains("printer") && json_inp["printer"].contains("file_name")) {
            base = detail::remove_extension(json_inp["printer"]["file_name"].get<std::string>());
        }
    } catch (...) {
        /* keep default */
    }

    if (mrcpp::mpi::grand_master()) {
        std::ofstream ofs(base + ".json", std::ios::out);
        ofs << json_tot.dump(2) << std::endl;
        ofs.close();
    }
}

// -----------------------------------------------------------------------------
// small helpers
// -----------------------------------------------------------------------------
std::string detail::remove_extension(const std::string &fname) {
    size_t lastdot = fname.find_last_of(".");
    if (lastdot == std::string::npos) return fname;
    return fname.substr(0, lastdot);
}

bool detail::all_success(const json &json_out) {
    auto scf_success = json_out["scf_calculation"]["success"].get<bool>();
    auto rsp_success = true;
    for (const auto &x : json_out["rsp_calculations"]) { rsp_success &= x["success"].get<bool>(); }
    return scf_success & rsp_success;
}

} // namespace mrchem
