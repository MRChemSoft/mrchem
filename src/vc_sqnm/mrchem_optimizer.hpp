#include <iostream>
#include <string>
#include<tuple>
#include <fstream>
#include <vector>

#include <MRCPP/Timer>
#include <MRCPP/Parallel>

#include "driver.h"
#include "mrchem.h"
#include "mrenv.h"
#include "version.h"

#include "chemistry/Molecule.h"
#include "chemistry/PhysicalConstants.h"
#include "vc_sqnm/periodic_optimizer.hpp"

#include <Eigen/Dense>
#include "MRCPP/Printer"

using json = nlohmann::json;
using namespace mrchem;

/**
 * @brief Writes positions to xyz file. Unit if file is angstrom.
 * @param xyzFile open file stream.
 * @param atomicPositions Eigen matrix containing the atomic positions, shape(3, n). Unit is bohr, will be converted to angstrom.
 * @param atomicLabels Vector containing strings of all the chemical symbols of all atoms
 * @param comment String containing comment that will be printed on the second line
*/
void writeXYZFile(std::ofstream& xyzFile, const Eigen::MatrixXd& atomicPositions, const std::vector<std::string>& atomicLabels, std::string comment) {
    if (!xyzFile.is_open()) {
        std::cerr << "Error: Output stream is not open for writing." << std::endl;
        return;
    }

    int numAtoms = atomicPositions.cols();

    double bohrToAng = 1.0 / PhysicalConstants::get("angstrom2bohrs");

    xyzFile << numAtoms << "\n";
    xyzFile << comment << std::endl;

    for (int i = 0; i < numAtoms; ++i) {
        xyzFile << atomicLabels[i] << " " << atomicPositions(0, i) * bohrToAng << " " << atomicPositions(1, i) * bohrToAng
                << " " << atomicPositions(2, i) * bohrToAng << "\n";
    }
    xyzFile.flush();
}

/**
 * @brief Gets postion from a molecule in json format.
 * @param mol_inp Molecule in json format.
 * 
 * @return Positions as an Eigen::matrix of shape (3, num_atoms)
*/
Eigen::MatrixXd getPositions(const json &mol_inp) {
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    Eigen::MatrixXd pos(3, mol.getNNuclei());

    int i = 0;
    for (const auto &coord : mol_inp["coords"]) {
        pos.col(i) << coord["xyz"][0], coord["xyz"][1], coord["xyz"][2];
        i++;
    }
    return pos;
}

/**
 * @brief Sets positions in json molecule.
 * @param mol_inp Inpunt molecule in json format.
 * @param pos Eigen::MatrixXd containing positions. Dimension must be (3, num_atoms)
*/
void setPositions(json &mol_inp, const Eigen::MatrixXd &pos) {
    for (int i = 0; i < mol_inp["coords"].size(); i++) {
        Eigen::VectorXd vvv = pos.col(i);
        mol_inp["coords"][i]["xyz"] = {vvv(0), vvv(1), vvv(2)};
    }
}

/**
 * @brief Does an scf calculation of a molecule.
 * 
 * @param mol_inp: json that describes the molecule.
 * @param scf_inp: scf settings.
 * 
 * @return: tuple containing print_properties of scf results the json that the driver returned
*/
std::tuple<json, json> getSCFResults(const json mol_inp, const json scf_inp) {
    Molecule mol;
    std::tuple <json, json> results;
    driver::init_molecule(mol_inp, mol);
    json scf_out = driver::scf::run(scf_inp, mol);
    // keeping the mpi barrier to be on the safe side, but not sure if it is needed
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    results = std::make_tuple(driver::print_properties(mol), scf_out);
    return results;
}

/**
 * @brief Extracts forces from json containing scf results.
 * @param scf_results json containing the results of an scf calculations. Must be obtained using the function getSCFResults.
 * @param forces: Provide empty matrix of dimension (3, num_atoms). Contains forces after function was called.
 * 
 * @return Forces acting on nuclei
*/
void extractForcesInPlace(const json &scf_results, Eigen::MatrixXd &forces){
    for (int i = 0; i < forces.cols(); i++)
    {
        for (int j = 0; j < 3; j++) {
            forces(j, i) = scf_results["geometric_derivative"]["geom-1"]["total"][3*i + j];
        }
    }
    // forces are the negative nuclear gradient.
    forces = - forces;
}

/**
 * @brief Extracts energy from json containing scf results.
 * @param scf_results json containing the results of an scf calculations. Must be obtained using the function getSCFResults.
 * 
 * @return Total energy
*/
double extractEnergy(const json &scf_results){
    return scf_results["scf_energy"]["E_tot"];
}

/**
 * @brief Optimizes positions of nuclei.
 * 
 * @param scf_inp: scf settings.
 * @param mol_inp: json that contains the molecule.
 * @param geopt_inp: json that contains the geometry optization settings.
 * @param jsonOutFileName the filename of the default json output. Used to construct names for own output files.
 * @return: A summary of the geometry optimization trajectory.
*/
json optimize_positions(json scf_inp, json mol_inp, const json &geopt_inp, std::string jsonOutFileName) {

    size_t stringPos = jsonOutFileName.find(".json");
    if (stringPos != std::string::npos)
    {
        jsonOutFileName.erase(stringPos, 5);
    }

    int num_atoms = mol_inp["coords"].size();
    int printLevel = 0;
    int pprec = 2 * mrcpp::Printer::getPrecision();

    std::ofstream xyzFile;
    std::vector<std::string> element_symbols;
    std::string element_symbol;

    for (int i = 0; i < mol_inp["coords"].size(); i++) {
        element_symbol = mol_inp["coords"][i]["atom"];
        element_symbol[0] = std::toupper(element_symbol[0]);
        element_symbols.push_back(element_symbol);
    }

    if (mrcpp::mpi::grand_master()) {
        xyzFile.open(jsonOutFileName + ".xyz", std::ios::out);
    }

    mrcpp::print::header(printLevel, "Starting geometry optimization using the SQNM method", 1, '=');
    println(printLevel, " Scientific users of the geometry optimization feature should cite");
    println(printLevel, " M. Gubler, M. Krummenacher, H. Huber, S. Goedecker");
    println(printLevel, " Journal of Computational Physics: X 2023, DOI: 10.1016/j.jcpx.2023.100131\n");
    mrcpp::print::separator(printLevel, '=', 2);

    // define default parameters
    // The sqnm parameters are documented in the periodic_optimizer.hpp file.
    auto max_iter = geopt_inp["max_iter"];
    auto max_history_length = geopt_inp["max_history_length"];
    auto init_step_size = geopt_inp["init_step_size"];
    auto minimal_step_size = geopt_inp["minimal_step_size"];
    auto subspace_tolerance = geopt_inp["subspace_tolerance"];
    auto max_force_component = geopt_inp["max_force_component"];

    PES_optimizer::periodic_optimizer optimizer(num_atoms, init_step_size, max_history_length, minimal_step_size, subspace_tolerance);
    
    int i = 0;
    Eigen::MatrixXd forces(3, num_atoms);
    double energy;
    
    std::tuple<json, json> results_tuple = getSCFResults(mol_inp, scf_inp);
    json results = std::get<0>(results_tuple);
    energy = extractEnergy(results);

    double energyOld = energy;
    extractForcesInPlace(results, forces);
    mrcpp::print::header(-1, "Geometry optimization summary of initial iteration", 0, '=');
    print_utils::scalar(0, "Iteration", i, "", 0);
    print_utils::scalar(-1, "Energy", energy, "Ha", pprec, true);
    print_utils::scalar(0, "Maximal force component", forces.cwiseAbs().maxCoeff(), "Ha / Bohr", pprec, true);
    print_utils::scalar(0, "Convergence threshold", max_force_component, "Ha / Bohr", pprec, true);
    mrcpp::print::separator(printLevel, '=', 2);

    json summary;
    summary["iteration_" + std::to_string(i)] = {
        {"results", results},
        {"molecule", mol_inp}
    };

    i++;

    Eigen::MatrixXd pos = getPositions(mol_inp);
    if (mrcpp::mpi::grand_master()) {
        std::string comment = "Unit of positions is angstrom. Iteration 0, total energy:" + std::to_string(energy) 
            + " max force component: " + std::to_string(forces.cwiseAbs().maxCoeff());
        writeXYZFile(xyzFile, pos, element_symbols, comment);
    }

    while (i < max_iter && forces.cwiseAbs().maxCoeff() > max_force_component) {
        optimizer.step(pos, energy, forces);
        setPositions(mol_inp, pos);
        if (geopt_inp["use_previous_guess"]) {
            scf_inp["initial_guess"]["type"] = "mw";
            scf_inp["initial_guess"]["file_phi_p"] = scf_inp["write_orbitals"]["file_phi_p"];
            scf_inp["initial_guess"]["file_phi_a"] = scf_inp["write_orbitals"]["file_phi_a"];
            scf_inp["initial_guess"]["file_phi_b"] = scf_inp["write_orbitals"]["file_phi_b"];
        }
        std::tuple<json, json> results_tuple = getSCFResults(mol_inp, scf_inp);
        json results = std::get<0>(results_tuple);
        energy = extractEnergy(results);
        extractForcesInPlace(results, forces);
        summary["iteration_" + std::to_string(i)] = {
            {"results", results},
            {"molecule", mol_inp},
            {"energy", energy},
            {"max_force_component", forces.cwiseAbs().maxCoeff()},
        };
        if (mrcpp::mpi::grand_master()) {
            std::string comment = "Iteration " + std::to_string(i) + ", total energy: " + std::to_string(energy)
                + " max force component: " + std::to_string(forces.cwiseAbs().maxCoeff());
            writeXYZFile(xyzFile, pos, element_symbols, comment);
        }
        mrcpp::print::header(printLevel, "Geometry optimization summary of iteration", 0, '=');
        print_utils::scalar(0, "Iteration", i, "", 0);
        print_utils::scalar(-1, "Energy", energy, "Ha", pprec, true);
        print_utils::scalar(0, "Energy improvement", energyOld - energy, "Ha", pprec, true);
        print_utils::scalar(0, "Maximal force component", forces.cwiseAbs().maxCoeff(), "Ha / Bohr", pprec, true);
        print_utils::scalar(0, "Convergence threshold", max_force_component, "Ha / Bohr", pprec, true);
        mrcpp::print::separator(printLevel, '=', 2);
        energyOld = energy;
        i++;
    }

    // print warning if geometry optimization did not converge
    if (i >= max_iter)
    {
        println(printLevel, "The geometry optimization did not converge!!!")
    }
    
    // make last step for correct ground state energy estimation.
    optimizer.step(pos, energy, forces);
    mrcpp::print::separator(-1, '=', 0);
    mrcpp::print::value(-1, "Est. energy of minimum", optimizer.lower_bound(), "Ha", pprec);
    mrcpp::print::value(-1, "Est. difference from minimum", energy - optimizer.lower_bound(), "Ha", pprec);
    mrcpp::print::separator(-1, '=', 0);

    if (mrcpp::mpi::grand_master()) {
        std::ofstream ofs;
        ofs.open(jsonOutFileName + "_optimization_summary.json", std::ios::out);
        ofs << summary.dump(2) << std::endl;
        ofs.close();
        xyzFile.close();
    }

    return std::get<1>(results_tuple);
}
