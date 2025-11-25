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

#include <iostream>
#include <string>

#include <MRCPP/Timer>
#include <MRCPP/Parallel>

#include "driver.h"
#include "mrchem.h"
#include "mrenv.h"
#include "version.h"

#include "chemistry/Molecule.h"
#include "chemistry/PhysicalConstants.h"

#include "vc_sqnm/mrchem_optimizer.hpp"

// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

int main(int argc, char **argv) {
    if (std::string(argv[1]) == "--version" || std::string(argv[1]) == "-v") {
        std::cout << program_version() << std::endl;
        return EXIT_SUCCESS;
    }

    const auto json_inp = mrenv::fetch_json(argc, argv);

    mrenv::initialize(json_inp);
    const auto &mol_inp = json_inp["molecule"];
    const auto &scf_inp = json_inp["scf_calculation"];
    const auto &rsp_inp = json_inp["rsp_calculations"];
    const auto &con_inp = json_inp["constants"];
    const auto &geopt_inp = json_inp["geom_opt"];

    // Instantiate the physical constants singleton
    PhysicalConstants::Initialize(con_inp);
    if (json_inp["printer"]["print_constants"]) PhysicalConstants::Print();

    Timer timer;
    json json_out;

    if (geopt_inp["run"]) {
        json_out = optimize_positions(scf_inp, mol_inp, geopt_inp, json_inp["printer"]["file_name"]);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    } else {
        Molecule mol;
        driver::init_molecule(mol_inp, mol);
        auto scf_out = driver::scf::run(scf_inp, mol);
        json rsp_out = {};
        if (scf_out["success"]) {
            for (auto &i : rsp_inp.items()) rsp_out[i.key()] = driver::rsp::run(i.value(), mol);
        }
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
        // Name and version of the output schema
        json_out["schema_name"] = "mrchem_output";
        json_out["schema_version"] = 1;
        // Computed values
        json_out["scf_calculation"] = scf_out;
        json_out["rsp_calculations"] = rsp_out;
        json_out["properties"] = driver::print_properties(mol);
        // Global success field: true if all requested calculations succeeded
        json_out["success"] = detail::all_success(json_out);
    }
    // Provenance
    json_out["provenance"] = {{"creator", "MRChem"},
                              {"version", program_version()},
                              {"nthreads", mrcpp::omp::n_threads},
                              {"mpi_processes", mrcpp::mpi::world_size},
                              {"total_cores", mrcpp::mpi::world_size * mrcpp::omp::n_threads},
                              {"routine", "mrchem.x"}};

    mrenv::finalize(timer.elapsed());
    mrenv::dump_json(json_inp, json_out);
    mrcpp::mpi::finalize();
    return EXIT_SUCCESS;
}
