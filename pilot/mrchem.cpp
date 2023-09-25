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

/** The MRChem sandbox */

#include "MRCPP/Parallel"
#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "chemistry/Nucleus.h"
#include "chemistry/PeriodicTable.h"
#include "environment/Cavity.h"
#include "environment/DHScreening.h"
#include "environment/GPESolver.h"
#include "environment/Permittivity.h"

#include "qmfunctions/qmfunction_fwd.h"
#include "qmoperators/two_electron/ReactionOperator.h"

#include "mrchem.h"
#include "mrenv.h"

// Initializing global variables
mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using mrcpp::ABGVOperator;

using json = nlohmann::json;
using Timer = mrcpp::Timer;
using namespace mrchem;

using DerivativeOperator = mrcpp::DerivativeOperator<3>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;

int main(int argc, char **argv) {
    mrcpp::mpi::initialize();
    const auto json_inp = mrenv::fetch_json(argc, argv);
    mrenv::initialize(json_inp);

    Timer timer;

    // Do your stuff here
    bool acc_pot = true;
    bool dyn_thrs = false;
    int kain = 7;
    int max_iter = 200;
    double poisson_prec = 1.0e-5;
    double eps_in = 1.0;
    double eps_out = 78.54;
    double kappa_out = 0.054995;
    double slope = 0.2;

    std::cout << "setting up the cavity" << std::endl;

    std::vector<double> R = {3.7794522509156563}; // must be 2 A
    auto sph_coords = std::vector<mrcpp::Coord<3>>({{0.0, 0.0, 0.0}});
    Cavity C(sph_coords, R, slope);

    std::cout << "setting up permittivity and kappa" << std::endl;
    Permittivity eps(C, eps_in, eps_out, "exponential");
    DHScreening kappa_sq(C, kappa_out, "continuous");

    std::cout << "setting up the operators" << std::endl;
    PoissonOperator_p P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
    DerivativeOperator_p D_p = std::make_shared<ABGVOperator<3>>(*MRA, 0.0, 0.0);

    std::cout << "setting up the nuclei" << std::endl;
    PeriodicTable P;
    auto q_coords = std::vector<mrcpp::Coord<3>>({{0.7558904498503081, 0.0, 0.0},
                                                  {0.0, 1.5117808997006161, 0.0},
                                                  {0.0, 0.0, 2.267671349550924},
                                                  {0.0, 0.0, -0.7558904498503081},
                                                  {-1.5117808997006161, 0.0, 0.0},
                                                  {0.0, -2.267671349550924, 0.0}});

    std::cout << "please work";
    Nucleus Q1(P.getElement(0), q_coords[0]);
    Nucleus Q2(P.getElement(0), q_coords[1]);
    Nucleus Q3(P.getElement(0), q_coords[2]);
    Nucleus Q4(P.getElement(0), q_coords[3]);
    Nucleus Q5(P.getElement(0), q_coords[4]);
    Nucleus Q6(P.getElement(0), q_coords[5]);

    Nuclei N;
    N.push_back(Q1);
    N.push_back(Q2);
    N.push_back(Q3);
    N.push_back(Q4);
    N.push_back(Q5);
    N.push_back(Q6);

    auto ORB = std::make_shared<mrchem::OrbitalVector>();
    ORB->push_back(Orbital(SPIN::Paired));

    std::cout << "setting up the solver" << std::endl;
    auto PB_solver = std::make_unique<GPESolver>(eps, N, P_p, D_p, poisson_prec, kain, max_iter, acc_pot, dyn_thrs, "nuclear");
    // auto PB_solver = std::make_unique<GPESolver>(eps, kappa_sq, N, P_p, D_p, poisson_prec, kain, max_iter, acc_pot, dyn_thrs, "nuclear");
    std::cout << "Computing the potential standard: " << std::endl;

    PB_solver->setConvergenceThreshold(poisson_prec);
    // PB_solver->setStaticSalt();
    PB_solver->iterateEquation(poisson_prec, ORB);
    double total_energy = PB_solver->getTotalEnergy() / 2;
    PB_solver->clear();

    total_energy = PB_solver->getTotalEnergy() / 2.0;
    std::cout << "total electrostatic potential linearized pb: " << total_energy << std::endl;

    println(0, json_inp.dump(2));

    timer.stop();
    double wt = timer.elapsed();

    mrenv::finalize(wt);
    mrcpp::mpi::finalize();
    return EXIT_SUCCESS;
}
