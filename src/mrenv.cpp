/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <MRCPP/Printer>
#include <XCFun/xcfun.h>
#include <fstream>

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"
#include "utils/print_utils.h"
#include "version.h"

using json = nlohmann::json;
using Printer = mrcpp::Printer;

namespace mrchem {

namespace mrenv {
void init_printer(const json &json_print);
void init_mra(const json &json_mra);
void init_mpi(const json &json_mpi);
void print_header();
} // namespace mrenv

json mrenv::fetch_json(int argc, char **argv) {
    const char *infile = nullptr;
    if (argc == 1) {
        infile = "STDIN";
    } else if (argc == 2) {
        infile = argv[1];
    } else {
        MSG_ERROR("Ivalid number of arguments!");
    }

    // Read JSON input
    json input;
    std::ifstream ifs(infile, std::ios_base::in);
    ifs >> input;
    ifs.close();

    return input["input"];
}

void mrenv::initialize(const json &json_inp) {
    auto json_print = json_inp.find("printer");
    auto json_mra = json_inp.find("mra");
    auto json_mpi = json_inp.find("mpi");

    if (json_mra == json_inp.end()) {
        MSG_ABORT("Missing MRA input!");
    } else {
        mrenv::init_mra(*json_mra);
    }
    if (json_mpi != json_inp.end()) mrenv::init_mpi(*json_mpi);
    if (json_print != json_inp.end()) mrenv::init_printer(*json_print);

    mrenv::print_header();
}

void mrenv::init_printer(const json &json_print) {
    // Initialize printing
    auto print_level = json_print["print_level"];
    auto print_prec = json_print["print_prec"];
    auto print_width = json_print["print_width"];
    auto print_mpi = json_print["print_mpi"];
    auto fname = json_print["file_name"].get<std::string>();
    if (print_mpi) {
        Printer::init(print_level, mpi::world_rank, mpi::world_size, fname.c_str());
    } else {
        Printer::init(print_level, mpi::world_rank, mpi::world_size);
    }
    Printer::setPrecision(print_prec);
    Printer::setWidth(print_width);
}

void mrenv::init_mra(const json &json_mra) {
    // Initialize world box
    int min_scale = json_mra["min_scale"];
    int max_scale = json_mra["max_scale"];
    auto corner = json_mra["corner"];
    auto boxes = json_mra["boxes"];
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes);

    // Initialize scaling basis
    auto order = json_mra["basis_order"];
    auto btype = json_mra["basis_type"];

    auto max_depth = max_scale - min_scale;
    if (min_scale < mrcpp::MinScale) MSG_ABORT("Root scale too large");
    if (max_scale > mrcpp::MaxScale) MSG_ABORT("Max scale too large");
    if (max_depth > mrcpp::MaxDepth) MSG_ABORT("Max depth too large");

    // Initialize global MRA
    if (btype == "interpolating") {
        mrcpp::InterpolatingBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else if (btype == "legendre") {
        mrcpp::LegendreBasis basis(order);
        MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
    } else {
        MSG_ABORT("Invalid basis type!");
    }
}

void mrenv::init_mpi(const json &json_mpi) {
    mpi::numerically_exact = json_mpi["numerically_exact"];
    mpi::shared_memory_size = json_mpi["shared_memory_size"];
    mpi::bank_size = json_mpi["bank_size"];
    mpi::initialize(); // NB: must be after bank_size and init_mra but before init_printer and print_header
}

void mrenv::print_header() {
    auto pwidth = Printer::getWidth();
    auto txt_width = 50;
    auto pre_spaces = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');
    std::stringstream o_ver, o_branch, o_hash, o_author, o_date;
    o_ver << "VERSION            " << program_version();
    o_branch << "Git branch         " << git_branch();
    o_hash << "Git commit hash    " << git_commit_hash();
    o_author << "Git commit author  " << git_commit_author();
    o_date << "Git commit date    " << git_commit_date();

    int ver_len = o_ver.str().size();
    int branch_len = o_branch.str().size();
    int hash_len = o_hash.str().size();
    int auth_len = o_author.str().size();
    int date_len = o_date.str().size();

    o_ver << std::string(std::max(0, txt_width - ver_len), ' ');
    o_branch << std::string(std::max(0, txt_width - branch_len), ' ');
    o_hash << std::string(std::max(0, txt_width - hash_len), ' ');
    o_author << std::string(std::max(0, txt_width - auth_len), ' ');
    o_date << std::string(std::max(0, txt_width - date_len), ' ');

    std::stringstream o_bank;
    if (mpi::bank_size > 0) {
        o_bank << "(" << mpi::bank_size << " bank)";
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
    print_utils::scalar(0, "MPI processes  ", mpi::world_size, o_bank.str(), 0, false);
    print_utils::scalar(0, "OpenMP threads ", omp::n_threads, "", 0, false);
    print_utils::scalar(0, "Total cores    ", mpi::world_size * omp::n_threads, "", 0, false);
    mrcpp::print::separator(0, ' ');
    mrcpp::print::separator(0, '-', 1);
    printout(0, xcfun_splash());
    mrcpp::print::environment(0);
    MRA->print();
}

void mrenv::finalize(double wt) {
    // Delete global MRA
    if (MRA != nullptr) delete MRA;
    MRA = nullptr;

    auto pwidth = Printer::getWidth();
    auto txt_width = 45;
    auto pre_spaces = (pwidth - 6 - txt_width) / 2;
    auto post_spaces = pwidth - 6 - txt_width - pre_spaces;
    std::string pre_str = std::string(3, '*') + std::string(pre_spaces, ' ');
    std::string post_str = std::string(post_spaces, ' ') + std::string(3, '*');

    auto hr = static_cast<int>(wt / 3600.0);
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
    json_tot["input"] = json_inp;
    json_tot["output"] = json_out;

    const auto file_name = detail::remove_extension(json_inp["printer"]["file_name"].get<std::string>());
    if (mpi::grand_master()) {
        std::ofstream ofs;
        ofs.open(file_name + ".json", std::ios::out);
        ofs << json_tot.dump(2) << std::endl;
        ofs.close();
    }
}

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
