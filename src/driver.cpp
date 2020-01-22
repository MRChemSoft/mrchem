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

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "driver.h"

#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"

#include "initial_guess/chk.h"
#include "initial_guess/core.h"
#include "initial_guess/gto.h"
#include "initial_guess/mw.h"
#include "initial_guess/sad.h"

#include "utils/MolPlotter.h"
#include "utils/math_utils.h"
#include "utils/print_utils.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/ExchangeOperator.h"
#include "qmoperators/two_electron/FockOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

#include "qmoperators/one_electron/H_BB_dia.h"
#include "qmoperators/one_electron/H_BM_dia.h"
#include "qmoperators/one_electron/H_B_dip.h"
#include "qmoperators/one_electron/H_B_spin.h"
#include "qmoperators/one_electron/H_E_dip.h"
#include "qmoperators/one_electron/H_E_quad.h"
#include "qmoperators/one_electron/H_M_fc.h"
#include "qmoperators/one_electron/H_M_pso.h"
#include "qmoperators/one_electron/NuclearGradientOperator.h"

#include "scf_solver/GroundStateSolver.h"
#include "scf_solver/KAIN.h"
#include "scf_solver/LinearResponseSolver.h"

#include "mrdft/Factory.h"

using mrcpp::ABGVOperator;
using mrcpp::Coord;
using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

using nlohmann::json;

using DerivativeOperator = mrcpp::DerivativeOperator<3>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

namespace mrchem {
namespace driver {

DerivativeOperator_p get_derivative(const std::string &name);
template <int I> RankOneTensorOperator<I> get_operator(const json &json_oper);
template <int I, int J> RankTwoTensorOperator<I, J> get_operator(const json &json_oper);
void build_fock_operator(const json &input, Molecule &mol, FockOperator &F, int order);

namespace scf {
bool guess_orbitals(const json &input, Molecule &mol);
bool guess_energy(const json &input, Molecule &mol, FockOperator &F);
void write_orbitals(const json &input, Molecule &mol);
void calc_properties(const json &input, Molecule &mol);
void plot_quantities(const json &input, Molecule &mol);
} // namespace scf

namespace rsp {
bool guess_orbitals(const json &input, Molecule &mol);
void write_orbitals(const json &input, Molecule &mol, bool dynamic);
void calc_properties(const json &input, Molecule &mol, int dir, double omega);
} // namespace rsp

} // namespace driver

/** @brief Initialize a molecule from input
 *
 * This function expects the "molecule" subsection of the input.
 */
void driver::init_molecule(const json &json_mol, Molecule &mol) {
    print_utils::headline(0, "Initializing Molecule");

    auto charge = json_mol["charge"];
    auto multiplicity = json_mol["multiplicity"];
    auto gauge_origin = json_mol["gauge_origin"];

    mol.setCharge(charge);
    mol.setMultiplicity(multiplicity);
    mol.setGaugeOrigin(gauge_origin);

    auto &nuclei = mol.getNuclei();
    for (const auto &coord : json_mol["coords"]) {
        auto atom = coord["atom"];
        auto xyz = coord["xyz"];
        nuclei.push_back(atom, xyz);
    }

    mol.printGeometry();
}

/** @brief Run ground-state SCF calculation
 *
 * This function will update the ground state orbitals and the Fock
 * matrix of the molecule based on the chosen electronic structure
 * method. The resulting orbitals will be either diagonalized or
 * localized at exit. Returns true if the calculation converges.
 *
 * After convergence the requested ground-state properties are computed.
 *
 * This function expects the "scf_calculation" subsection of the input.
 */
json driver::scf::run(const json &json_scf, Molecule &mol) {
    print_utils::headline(0, "Computing Ground State Wavefunction");
    json json_out;
    json_out["success"] = true;

    ///////////////////////////////////////////////////////////
    ////////////////   Building Fock Operator   ///////////////
    ///////////////////////////////////////////////////////////

    FockOperator F;
    const auto &json_fock = json_scf["fock_operator"];
    driver::build_fock_operator(json_fock, mol, F, 0);

    ///////////////////////////////////////////////////////////
    ///////////////   Setting Up Initial Guess   //////////////
    ///////////////////////////////////////////////////////////

    const auto &json_guess = json_scf["initial_guess"];
    if (driver::scf::guess_orbitals(json_guess, mol)) {
        driver::scf::guess_energy(json_guess, mol, F);
        json_out["initial_energy"] = mol.getSCFEnergy().json();
    } else {
        json_out["success"] = false;
        return json_out;
    }

    ///////////////////////////////////////////////////////////
    //////////   Optimizing Ground State Orbitals  ////////////
    ///////////////////////////////////////////////////////////

    // Run GroundStateSolver if present in input JSON
    auto scf_solver = json_scf.find("scf_solver");
    if (scf_solver != json_scf.end()) {
        auto kain = (*scf_solver)["kain"];
        auto method = (*scf_solver)["method"];
        auto max_iter = (*scf_solver)["max_iter"];
        auto rotation = (*scf_solver)["rotation"];
        auto localize = (*scf_solver)["localize"];
        auto file_chk = (*scf_solver)["file_chk"];
        auto checkpoint = (*scf_solver)["checkpoint"];
        auto start_prec = (*scf_solver)["start_prec"];
        auto final_prec = (*scf_solver)["final_prec"];
        auto energy_thrs = (*scf_solver)["energy_thrs"];
        auto orbital_thrs = (*scf_solver)["orbital_thrs"];
        auto helmholtz_prec = (*scf_solver)["helmholtz_prec"];

        GroundStateSolver solver;
        solver.setHistory(kain);
        solver.setRotation(rotation);
        solver.setLocalize(localize);
        solver.setMethodName(method);
        solver.setCheckpoint(checkpoint);
        solver.setCheckpointFile(file_chk);
        solver.setMaxIterations(max_iter);
        solver.setHelmholtzPrec(helmholtz_prec);
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, energy_thrs);
        json_out["scf_solver"] = solver.optimize(mol, F);
        json_out["success"] = json_out["scf_solver"]["converged"];
    }

    ///////////////////////////////////////////////////////////
    //////////   Computing Ground State Properties   //////////
    ///////////////////////////////////////////////////////////

    if (json_out["success"]) {
        auto json_orbs = json_scf.find("write_orbitals");
        if (json_orbs != json_scf.end()) driver::scf::write_orbitals(*json_orbs, mol);

        auto json_prop = json_scf.find("properties");
        if (json_prop != json_scf.end()) driver::scf::calc_properties(*json_prop, mol);

        auto json_plot = json_scf.find("cube_plot");
        if (json_plot != json_scf.end()) driver::scf::plot_quantities(*json_plot, mol);
    }

    return json_out;
}

/** @brief Run initial guess calculation for the orbitals
 *
 * This function will update the ground state orbitals and the Fock
 * matrix of the molecule, based on the chosen initial guess method.
 * The orbital vector is initialized with the appropriate particle
 * number and spin. The Fock matrix is initialized to the zero matrix
 * of the appropriate size.
 *
 * This function expects the "initial_guess" subsection of the input.
 */
bool driver::scf::guess_orbitals(const json &json_guess, Molecule &mol) {
    auto prec = json_guess["prec"];
    auto zeta = json_guess["zeta"];
    auto type = json_guess["type"];
    auto mw_p = json_guess["file_phi_p"];
    auto mw_a = json_guess["file_phi_a"];
    auto mw_b = json_guess["file_phi_b"];
    auto gto_p = json_guess["file_gto_p"];
    auto gto_a = json_guess["file_gto_a"];
    auto gto_b = json_guess["file_gto_b"];
    auto gto_bas = json_guess["file_basis"];
    auto file_chk = json_guess["file_chk"];
    auto restricted = json_guess["restricted"];

    // Figure out number of electrons
    int mult = mol.getMultiplicity(); // multiplicity
    int Ne = mol.getNElectrons();     // total electrons
    int Ns = mult - 1;                // single occ electrons
    int Nd = Ne - Ns;                 // double occ electrons
    if (Nd % 2 != 0) {
        MSG_ERROR("Invalid multiplicity");
        return false;
    }

    // Figure out number of occupied orbitals
    int Na = (restricted) ? Ns : Nd / 2 + Ns; // alpha orbitals
    int Nb = (restricted) ? 0 : Nd / 2;       // beta orbitals
    int Np = (restricted) ? Nd / 2 : 0;       // paired orbitals

    // Fill orbital vector
    auto &nucs = mol.getNuclei();
    auto &Phi = mol.getOrbitals();
    for (auto p = 0; p < Np; p++) Phi.push_back(Orbital(SPIN::Paired));
    for (auto a = 0; a < Na; a++) Phi.push_back(Orbital(SPIN::Alpha));
    for (auto b = 0; b < Nb; b++) Phi.push_back(Orbital(SPIN::Beta));
    mpi::distribute(Phi);

    auto success = true;
    if (type == "chk") {
        success = initial_guess::chk::setup(Phi, file_chk);
    } else if (type == "mw") {
        success = initial_guess::mw::setup(Phi, prec, mw_p, mw_a, mw_b);
    } else if (type == "core") {
        success = initial_guess::core::setup(Phi, prec, nucs, zeta);
    } else if (type == "sad") {
        success = initial_guess::sad::setup(Phi, prec, nucs, zeta);
    } else if (type == "gto") {
        success = initial_guess::gto::setup(Phi, prec, gto_bas, gto_p, gto_a, gto_b);
    } else {
        MSG_ERROR("Invalid initial guess");
        success = false;
    }
    orbital::print(Phi);

    return success;
}

bool driver::scf::guess_energy(const json &json_guess, Molecule &mol, FockOperator &F) {
    double prec = json_guess["prec"];
    auto method = json_guess["method"];
    auto localize = json_guess["localize"];

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation  ", "Compute initial energy");
    print_utils::text(0, "Method       ", method);
    print_utils::text(0, "Precision    ", print_utils::dbl_to_str(prec, 5, true));
    print_utils::text(0, "Localization ", (localize) ? "On" : "Off");
    mrcpp::print::separator(0, '~', 2);

    Timer timer;
    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "Calculating Molecular Energy");

    auto &Phi = mol.getOrbitals();
    auto &nucs = mol.getNuclei();
    auto &F_mat = mol.getFockMatrix();
    F_mat = ComplexMatrix::Zero(Phi.size(), Phi.size());
    if (localize) orbital::localize(prec, Phi, F_mat);

    F.setup(prec);
    F_mat = F(Phi, Phi);
    mol.getSCFEnergy() = F.trace(Phi, nucs);
    F.clear();
    if (plevel == 1) mrcpp::print::footer(1, timer, 2);
    if (not localize) orbital::diagonalize(prec, Phi, F_mat);
    mol.getSCFEnergy().print();

    return true;
}

void driver::scf::write_orbitals(const json &json_orbs, Molecule &mol) {
    auto &Phi = mol.getOrbitals();
    orbital::save_orbitals(Phi, json_orbs["file_phi_p"], SPIN::Paired);
    orbital::save_orbitals(Phi, json_orbs["file_phi_a"], SPIN::Alpha);
    orbital::save_orbitals(Phi, json_orbs["file_phi_b"], SPIN::Beta);
}

/** @brief Compute ground-state properties
 *
 * This function expects the "properties" subsection of the "scf_calculation"
 * input section, and will compute all properties which are present in this input.
 * This includes the diamagnetic contributions to the magnetic response properties.
 */
void driver::scf::calc_properties(const json &json_prop, Molecule &mol) {
    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "Computing Ground State Properties");

    auto &nuclei = mol.getNuclei();
    auto &Phi = mol.getOrbitals();
    auto &F_mat = mol.getFockMatrix();

    auto json_eps = json_prop.find("orbital_energies");
    if (json_eps != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing orbital energies");
        OrbitalEnergies &eps = mol.getOrbitalEnergies();
        eps.getOccupancy() = orbital::get_occupancies(Phi);
        eps.getEpsilon() = orbital::calc_eigenvalues(Phi, F_mat);
        eps.getSpin() = orbital::get_spins(Phi);
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Orbital energies", t_lap);
    }

    auto json_dip = json_prop.find("dipole_moment");
    if (json_dip != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing dipole moment");
        auto prec = (*json_dip)["precision"].get<double>();
        DipoleMoment &mu = mol.getDipoleMoment();

        auto h = driver::get_operator<3>(*json_dip);
        h.setup(prec);
        mu.getNuclear() = -h.trace(nuclei).real();
        mu.getElectronic() = h.trace(Phi).real();
        h.clear();
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Dipole moment", t_lap);
    }

    auto json_quad = json_prop.find("quadrupole_moment");
    if (json_quad != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing quadrupole moment");
        auto prec = (*json_quad)["precision"].get<double>();
        QuadrupoleMoment &Q = mol.getQuadrupoleMoment();

        auto h = driver::get_operator<3, 3>(*json_quad);
        h.setup(prec);
        Q.getNuclear() = -h.trace(nuclei).real();
        Q.getElectronic() = h.trace(Phi).real();
        h.clear();
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Quadrupole moment", t_lap);
    }

    auto json_hyperpolarizability = json_prop.find("hyperpolarizability");
    if (json_hyperpolarizability != json_prop.end()) MSG_ERROR("Hyperpolarizability not implemented");

    auto json_gradient = json_prop.find("nuclear_gradient");
    if (json_gradient != json_prop.end()) MSG_ERROR("Nuclear gradient not implemented");

    auto json_mag = json_prop.find("magnetizability");
    if (json_mag != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing magnetizability (dia)");
        auto prec = (*json_mag)["precision"].get<double>();
        Magnetizability &xi = mol.getMagnetizability();

        auto h = driver::get_operator<3, 3>(*json_mag);
        h.setup(prec);
        xi.getDiamagnetic() = -h.trace(Phi).real();
        h.clear();
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Magnetizability (dia)", t_lap);
    }

    auto json_nmr = json_prop.find("nmr_shielding");
    if (json_nmr != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing NMR shielding (dia)");
        for (const auto &json_nuc : *json_nmr) {
            auto k = json_nuc["nucleus_k"].get<int>();
            auto prec = json_nuc["precision"].get<double>();
            NMRShielding &sigma_k = mol.getNMRShielding(k);

            auto h = driver::get_operator<3, 3>(json_nuc);
            h.setup(prec);
            sigma_k.getDiamagnetic() = h.trace(Phi).real();
            h.clear();
        }
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "NMR shielding (dia)", t_lap);
    }
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
}

/** @brief Plot ground-state quantities
 *
 * This function expects the "cube_plot" subsection of the
 * "scf_calculation" input section.
 */
void driver::scf::plot_quantities(const json &json_plot, Molecule &mol) {
    Timer t_tot, t_lap;

    auto path = json_plot["plotter"]["path_plots"].get<std::string>();
    auto npts = json_plot["plotter"]["points"];
    auto O = json_plot["plotter"]["O"];
    auto A = json_plot["plotter"]["A"];
    auto B = json_plot["plotter"]["B"];
    auto C = json_plot["plotter"]["C"];
    auto dens_plot = json_plot["density"];
    auto orb_idx = json_plot["orbital"];

    print_utils::headline(1, "Plotting Ground State Quantities");
    mrcpp::print::header(1, "CubePlot");

    auto &Phi = mol.getOrbitals();
    MolPlotter plt(mol, O);
    plt.setRange(A, B, C);

    if (dens_plot) {
        Density rho(false);

        t_lap.start();
        std::string fname = path + "/rho_t";
        density::compute(-1.0, rho, Phi, DensityType::Total);
        plt.cubePlot(npts, rho, fname);
        rho.free(NUMBER::Total);
        mrcpp::print::time(1, fname, t_lap);

        if (orbital::size_singly(Phi) > 0) {
            t_lap.start();
            fname = path + "/rho_s";
            density::compute(-1.0, rho, Phi, DensityType::Spin);
            plt.cubePlot(npts, rho, fname);
            mrcpp::print::time(1, fname, t_lap);
            rho.free(NUMBER::Total);

            t_lap.start();
            fname = path + "/rho_a";
            density::compute(-1.0, rho, Phi, DensityType::Alpha);
            plt.cubePlot(npts, rho, fname);
            mrcpp::print::time(1, fname, t_lap);
            rho.free(NUMBER::Total);

            t_lap.start();
            fname = path + "/rho_b";
            density::compute(-1.0, rho, Phi, DensityType::Beta);
            plt.cubePlot(npts, rho, fname);
            rho.free(NUMBER::Total);
            mrcpp::print::time(1, fname, t_lap);
        }
    }

    // Plotting NO orbitals
    if (orb_idx.size() > 0) {
        if (orb_idx[0] < 0) {
            // Plotting ALL orbitals
            for (auto i = 0; i < Phi.size(); i++) {
                if (not mpi::my_orb(Phi[i])) continue;
                t_lap.start();
                std::stringstream name;
                name << path << "/phi_" << i;
                plt.cubePlot(npts, Phi[i], name.str());
                mrcpp::print::time(1, name.str(), t_lap);
            }
        } else {
            // Plotting some orbitals
            for (auto &i : orb_idx) {
                if (not mpi::my_orb(Phi[i])) continue;
                t_lap.start();
                std::stringstream name;
                name << path << "/phi_" << i;
                plt.cubePlot(npts, Phi[i], name.str());
                mrcpp::print::time(1, name.str(), t_lap);
            }
        }
    }

    mrcpp::print::footer(1, t_tot, 2);
}

/** @brief Run linear response SCF calculation
 *
 * This function will update the perturbed orbitals of the molecule
 * based on the chosen electronic structure method and perturbation
 * operator. Each response calculation corresponds to one particular
 * perturbation operator (could be a vector operator with several
 * components). Returns true if the calculation converges.
 *
 * After convergence the requested linear response properties are computed.
 *
 * This function expects a single subsection entry in the "rsp_calculations"
 * vector of the input.
 */
json driver::rsp::run(const json &json_rsp, Molecule &mol) {
    print_utils::headline(0, "Computing Linear Response Wavefunction");
    json json_out = {};

    ///////////////////////////////////////////////////////////
    /////////////   Preparing Unperturbed System   ////////////
    ///////////////////////////////////////////////////////////

    const auto &json_unpert = json_rsp["unperturbed"];
    const auto &unpert_fock = json_unpert["fock_operator"];
    auto unpert_loc = json_unpert["localize"];
    auto unpert_prec = json_unpert["prec"];

    auto &Phi = mol.getOrbitals();
    auto &F_mat = mol.getFockMatrix();

    if (unpert_loc) {
        orbital::localize(unpert_prec, Phi, F_mat);
    } else {
        orbital::diagonalize(unpert_prec, Phi, F_mat);
    }

    FockOperator F_0;
    driver::build_fock_operator(unpert_fock, mol, F_0, 0);
    F_0.setup(unpert_prec);

    ///////////////////////////////////////////////////////////
    //////////////   Preparing Perturbed System   /////////////
    ///////////////////////////////////////////////////////////

    auto omega = json_rsp["frequency"];
    auto dynamic = json_rsp["dynamic"];
    mol.initPerturbedOrbitals(dynamic);

    FockOperator F_1;
    const auto &json_fock_1 = json_rsp["fock_operator"];
    driver::build_fock_operator(json_fock_1, mol, F_1, 1);

    const auto &json_pert = json_rsp["perturbation"];
    auto h_1 = driver::get_operator<3>(json_pert);

    auto success = true;
    for (auto d = 0; d < 3; d++) {
        const auto &json_comp = json_rsp["components"][d];
        F_1.perturbation() = h_1[d];

        ///////////////////////////////////////////////////////////
        ///////////////   Setting Up Initial Guess   //////////////
        ///////////////////////////////////////////////////////////

        const auto &json_guess = json_comp["initial_guess"];
        success = driver::rsp::guess_orbitals(json_guess, mol);

        ///////////////////////////////////////////////////////////
        /////////////   Optimizing Perturbed Orbitals  ////////////
        ///////////////////////////////////////////////////////////

        auto rsp_solver = json_comp.find("rsp_solver");
        if (rsp_solver != json_comp.end()) {
            auto kain = (*rsp_solver)["kain"];
            auto method = (*rsp_solver)["method"];
            auto max_iter = (*rsp_solver)["max_iter"];
            auto file_chk_x = (*rsp_solver)["file_chk_x"];
            auto file_chk_y = (*rsp_solver)["file_chk_y"];
            auto checkpoint = (*rsp_solver)["checkpoint"];
            auto orth_prec = (*rsp_solver)["orth_prec"];
            auto start_prec = (*rsp_solver)["start_prec"];
            auto final_prec = (*rsp_solver)["final_prec"];
            auto orbital_thrs = (*rsp_solver)["orbital_thrs"];
            auto property_thrs = (*rsp_solver)["property_thrs"];
            auto helmholtz_prec = (*rsp_solver)["helmholtz_prec"];

            LinearResponseSolver solver(dynamic);
            solver.setHistory(kain);
            solver.setMethodName(method);
            solver.setMaxIterations(max_iter);
            solver.setCheckpoint(checkpoint);
            solver.setCheckpointFile(file_chk_x, file_chk_y);
            solver.setHelmholtzPrec(helmholtz_prec);
            solver.setOrbitalPrec(start_prec, final_prec);
            solver.setThreshold(orbital_thrs, property_thrs);
            solver.setOrthPrec(orth_prec);

            success = solver.optimize(omega, mol, F_0, F_1);
        }

        ///////////////////////////////////////////////////////////
        ////////////   Compute Response Properties   //////////////
        ///////////////////////////////////////////////////////////

        if (success) {
            auto json_orbs = json_comp.find("write_orbitals");
            if (json_orbs != json_comp.end()) driver::rsp::write_orbitals(*json_orbs, mol, dynamic);

            auto json_prop = json_rsp.find("properties");
            if (json_prop != json_rsp.end()) driver::rsp::calc_properties(*json_prop, mol, d, omega);
        }
        mol.getOrbitalsX().clear(); // Clear orbital vector
        mol.getOrbitalsY().clear(); // Clear orbital vector
    }
    F_0.clear();
    mol.getOrbitalsX_p().reset(); // Release shared_ptr
    mol.getOrbitalsY_p().reset(); // Release shared_ptr

    return json_out;
}

/** @brief Run initial guess calculation for the response orbitals
 *
 * This function will update the ground state orbitals and the Fock
 * matrix of the molecule, based on the chosen initial guess method.
 * The orbital vector is initialized with the appropriate particle
 * number and spin. The Fock matrix is initialized to the zero matrix
 * of the appropriate size.
 *
 * This function expects the "initial_guess" subsection of the input.
 */
bool driver::rsp::guess_orbitals(const json &json_guess, Molecule &mol) {
    auto prec = json_guess["prec"];
    auto type = json_guess["type"];
    auto mw_xp = json_guess["file_x_p"];
    auto mw_xa = json_guess["file_x_a"];
    auto mw_xb = json_guess["file_x_b"];
    auto mw_yp = json_guess["file_y_p"];
    auto mw_ya = json_guess["file_y_a"];
    auto mw_yb = json_guess["file_y_b"];
    auto file_chk_x = json_guess["file_chk_x"];
    auto file_chk_y = json_guess["file_chk_y"];

    auto &Phi = mol.getOrbitals();
    auto &X = mol.getOrbitalsX();
    auto &Y = mol.getOrbitalsY();

    auto success_x = false;
    X = orbital::param_copy(Phi);
    if (type == "chk") {
        success_x = initial_guess::chk::setup(X, file_chk_x);
    } else if (type == "mw") {
        success_x = initial_guess::mw::setup(X, prec, mw_xp, mw_xa, mw_xb);
    } else if (type == "none") {
        mrcpp::print::separator(0, '~');
        print_utils::text(0, "Calculation     ", "Compute initial orbitals");
        print_utils::text(0, "Method          ", "Zero guess");
        mrcpp::print::separator(0, '~', 2);
    } else {
        MSG_ERROR("Invalid initial guess");
    }
    orbital::print(X);

    auto success_y = false;
    if (&X != &Y) {
        Y = orbital::param_copy(Phi);
        if (type == "chk") {
            success_y = initial_guess::chk::setup(Y, file_chk_y);
        } else if (type == "mw") {
            success_y = initial_guess::mw::setup(Y, prec, mw_yp, mw_ya, mw_yb);
        } else if (type == "none") {
            mrcpp::print::separator(0, '~');
            print_utils::text(0, "Calculation     ", "Compute initial orbitals");
            print_utils::text(0, "Method          ", "Zero guess");
            mrcpp::print::separator(0, '~', 2);
        } else {
            MSG_ERROR("Invalid initial guess");
        }
        orbital::print(Y);
    }

    return (success_x and success_y);
}

void driver::rsp::write_orbitals(const json &json_orbs, Molecule &mol, bool dynamic) {
    auto &X = mol.getOrbitalsX();
    orbital::save_orbitals(X, json_orbs["file_x_p"], SPIN::Paired);
    orbital::save_orbitals(X, json_orbs["file_x_a"], SPIN::Alpha);
    orbital::save_orbitals(X, json_orbs["file_x_b"], SPIN::Beta);
    if (dynamic) {
        auto &Y = mol.getOrbitalsY();
        orbital::save_orbitals(Y, json_orbs["file_y_p"], SPIN::Paired);
        orbital::save_orbitals(Y, json_orbs["file_y_a"], SPIN::Alpha);
        orbital::save_orbitals(Y, json_orbs["file_y_b"], SPIN::Beta);
    }
}

/** @brief Compute linear response properties
 *
 * This function expects the "properties" subsection of the "rsp_calculations"
 * input section, and will compute all properties which are present in this input.
 */
void driver::rsp::calc_properties(const json &json_prop, Molecule &mol, int dir, double omega) {
    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "Computing Linear Response Properties");

    auto &Phi = mol.getOrbitals();
    auto &X = mol.getOrbitalsX();
    auto &Y = mol.getOrbitalsY();

    auto json_pol = json_prop.find("polarizability");
    if (json_pol != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing polarizability");
        auto prec = (*json_pol)["precision"].get<double>();
        Polarizability &alpha = mol.getPolarizability(omega);

        auto h = driver::get_operator<3>(*json_pol);
        h.setup(prec);
        alpha.getTensor().row(dir) = -h.trace(Phi, X, Y).real();
        h.clear();
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Polarizability", t_lap);
    }

    auto json_mag = json_prop.find("magnetizability");
    if (json_mag != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing magnetizability (para)");
        auto prec = (*json_mag)["precision"].get<double>();
        Magnetizability &xi = mol.getMagnetizability();

        auto h = driver::get_operator<3>(*json_mag);
        h.setup(prec);
        xi.getParamagnetic().row(dir) = -h.trace(Phi, X, Y).real();
        h.clear();
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Magnetizability (para)", t_lap);
    }

    auto json_nmr = json_prop.find("nmr_shielding");
    if (json_nmr != json_prop.end()) {
        t_lap.start();
        mrcpp::print::header(2, "Computing NMR shielding (para)");
        for (const auto &json_nuc : *json_nmr) {
            auto k = json_nuc["nucleus_k"].get<int>();
            auto prec = json_nuc["precision"].get<double>();
            NMRShielding &sigma_k = mol.getNMRShielding(k);

            auto h = driver::get_operator<3>(json_nuc);
            h.setup(prec);
            sigma_k.getParamagnetic().row(dir) = -h.trace(Phi, X, Y).real();
            h.clear();
        }
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "NMR shielding (para)", t_lap);
    }
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
}

/** @brief Build Fock operator based on input parameters
 *
 * This function expects the "fock_operator" subsection of input, and will
 * construct all operator which are present in this input. Option to set
 * perturbation order of the operators.
 */
void driver::build_fock_operator(const json &json_fock, Molecule &mol, FockOperator &F, int order) {
    auto &nuclei = mol.getNuclei();
    auto Phi_p = mol.getOrbitals_p();
    auto X_p = mol.getOrbitalsX_p();
    auto Y_p = mol.getOrbitalsY_p();

    ///////////////////////////////////////////////////////////
    //////////////////   Kinetic Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_kinetic = json_fock.find("kinetic_operator");
    if (json_kinetic != json_fock.end()) {
        auto kin_diff = (*json_kinetic)["derivative"].get<std::string>();
        auto D_p = driver::get_derivative(kin_diff);
        auto T_p = std::make_shared<KineticOperator>(D_p);
        F.getKineticOperator() = T_p;
    }
    ///////////////////////////////////////////////////////////
    //////////////////   Nuclear Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_nuclear = json_fock.find("nuclear_operator");
    if (json_nuclear != json_fock.end()) {
        auto proj_prec = (*json_nuclear)["proj_prec"].get<double>();
        auto smooth_prec = (*json_nuclear)["smooth_prec"].get<double>();
        auto shared_memory = (*json_nuclear)["shared_memory"].get<bool>();
        auto V_p = std::make_shared<NuclearOperator>(nuclei, proj_prec, smooth_prec, shared_memory);
        F.getNuclearOperator() = V_p;
    }
    ///////////////////////////////////////////////////////////
    //////////////////   Coulomb Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_coulomb = json_fock.find("coulomb_operator");
    if (json_coulomb != json_fock.end()) {
        auto poisson_prec = (*json_coulomb)["poisson_prec"].get<double>();
        auto shared_memory = (*json_coulomb)["shared_memory"].get<bool>();
        auto P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
        if (order == 0) {
            auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, shared_memory);
            F.getCoulombOperator() = J_p;
        } else if (order == 1) {
            auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, X_p, Y_p, shared_memory);
            F.getCoulombOperator() = J_p;
        } else {
            MSG_ABORT("Invalid perturbation order");
        }
    }
    ///////////////////////////////////////////////////////////
    ////////////////////   XC Operator   //////////////////////
    ///////////////////////////////////////////////////////////
    double exx = 1.0;
    auto json_xc = json_fock.find("xc_operator");
    if (json_xc != json_fock.end()) {
        auto grid_prec = (*json_xc)["grid_prec"].get<double>();
        auto shared_memory = (*json_xc)["shared_memory"].get<bool>();
        auto json_xcfunc = (*json_xc)["xc_functional"].get<json>();
        auto xc_spin = json_xcfunc["spin"].get<bool>();
        auto xc_gamma = json_xcfunc["gamma"].get<bool>();
        auto xc_log_grad = json_xcfunc["log_grad"].get<bool>();
        auto xc_cutoff = json_xcfunc["cutoff"].get<double>();
        auto xc_diff = json_xcfunc["derivative"].get<std::string>();
        auto xc_funcs = json_xcfunc["functionals"].get<json>();
        auto xc_order = order + 1;

        mrdft::Factory xc_factory(*MRA);
        xc_factory.setSpin(xc_spin);
        xc_factory.setOrder(xc_order);
        xc_factory.setUseGamma(xc_gamma);
        xc_factory.setLogGradient(xc_log_grad);
        xc_factory.setDensityCutoff(xc_cutoff);
        for (const auto &f : xc_funcs) {
            auto name = f["name"].get<std::string>();
            auto coef = f["coef"].get<double>();
            xc_factory.setFunctional(name, coef);
        }
        auto mrdft_p = xc_factory.build();
        exx = mrdft_p->functional().amountEXX();

        if (order == 0) {
            auto XC_p = std::make_shared<XCOperator>(mrdft_p, Phi_p, shared_memory);
            F.getXCOperator() = XC_p;
        } else if (order == 1) {
            auto XC_p = std::make_shared<XCOperator>(mrdft_p, Phi_p, X_p, Y_p, shared_memory);
            F.getXCOperator() = XC_p;
        } else {
            MSG_ABORT("Invalid perturbation order");
        }
    }
    ///////////////////////////////////////////////////////////
    /////////////////   Exchange Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_exchange = json_fock.find("exchange_operator");
    if (json_exchange != json_fock.end() and exx > mrcpp::MachineZero) {
        auto poisson_prec = (*json_exchange)["poisson_prec"].get<double>();
        auto screen_prec = (*json_exchange)["screen"].get<bool>();
        auto P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
        if (order == 0) {
            auto K_p = std::make_shared<ExchangeOperator>(P_p, Phi_p, screen_prec);
            F.getExchangeOperator() = K_p;
        } else {
            auto K_p = std::make_shared<ExchangeOperator>(P_p, Phi_p, X_p, Y_p, screen_prec);
            F.getExchangeOperator() = K_p;
        }
    }
    ///////////////////////////////////////////////////////////
    /////////////////   External Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_external = json_fock.find("external_operator");
    if (json_external != json_fock.end()) {
        auto field = (*json_external)["electric_field"].get<std::array<double, 3>>();
        auto r_O = (*json_external)["r_O"].get<Coord<3>>();
        auto V_ext = std::make_shared<ElectricFieldOperator>(field, r_O);
        F.getExtOperator() = V_ext;
    }
    F.build(exx);
}

/** @brief Construct perturbation operator based on input keyword */
template <int I> RankOneTensorOperator<I> driver::get_operator(const json &json_oper) {
    RankOneTensorOperator<I> h;
    auto oper = json_oper["operator"].get<std::string>();
    if (oper == "h_e_dip") {
        auto r_O = json_oper["r_O"].get<mrcpp::Coord<3>>();
        h = H_E_dip(r_O);
    } else if (oper == "h_b_dip") {
        auto r_O = json_oper["r_O"].get<mrcpp::Coord<3>>();
        auto pert_diff = json_oper["derivative"].get<std::string>();
        auto D = driver::get_derivative(pert_diff);
        h = H_B_dip(D, r_O);
    } else if (oper == "h_m_pso") {
        auto r_K = json_oper["r_K"].get<mrcpp::Coord<3>>();
        auto smoothing = json_oper["smoothing"].get<double>();
        auto pert_diff = json_oper["derivative"].get<std::string>();
        auto D = driver::get_derivative(pert_diff);
        h = H_M_pso(D, r_K, smoothing);
    } else {
        MSG_ERROR("Invalid operator: " << oper);
    }
    return h;
}

template <int I, int J> RankTwoTensorOperator<I, J> driver::get_operator(const json &json_oper) {
    RankTwoTensorOperator<I, J> h;
    auto oper = json_oper["operator"].get<std::string>();
    if (oper == "h_e_quad") {
        auto r_O = json_oper["r_O"].get<mrcpp::Coord<3>>();
        h = H_E_quad(r_O);
    } else if (oper == "h_bb_dia") {
        auto r_O = json_oper["r_O"].get<mrcpp::Coord<3>>();
        h = H_BB_dia(r_O);
    } else if (oper == "h_bm_dia") {
        auto r_O = json_oper["r_O"].get<mrcpp::Coord<3>>();
        auto r_K = json_oper["r_K"].get<mrcpp::Coord<3>>();
        auto smoothing = json_oper["smoothing"].get<double>();
        h = H_BM_dia(r_O, r_K, smoothing);
    } else {
        MSG_ERROR("Invalid operator: " << oper);
    }
    return h;
}

/** @brief Construct derivative operator based on input keyword */
DerivativeOperator_p driver::get_derivative(const std::string &name) {
    DerivativeOperator_p D = nullptr;
    if (name == "abgv_00") {
        D = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    } else if (name == "abgv_55") {
        D = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
    } else if (name == "ph") {
        D = std::make_shared<mrcpp::PHOperator<3>>(*MRA, 1);
    } else if (name == "bspline") {
        D = std::make_shared<mrcpp::BSOperator<3>>(*MRA, 1);
    } else {
        MSG_ERROR("Invalid derivative operator");
    }
    return D;
}

json driver::print_properties(const Molecule &mol) {
    print_utils::headline(0, "Printing Molecular Properties");
    mol.printGeometry();
    mol.printProperties();
    return mol.json();
}

} // namespace mrchem
