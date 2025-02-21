
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

#include <MRCPP/MWOperators>
#include <MRCPP/MWFunctions>
#include <MRCPP/Parallel>
#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <MRCPP/utils/details.h>

#include "core.h"
#include "gto.h"
#include "nao.h"

#include "chemistry/Nucleus.h"
#include "chemistry/chemistry_utils.h"
#include "pseudopotential/projectorOperator.h"
#include "pseudopotential/sphericalHarmonics.h"
#include "utils/PolyInterpolator.h"

#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

#include "utils/print_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/XCOperator.h"
#include "qmoperators/qmoperator_utils.h"

#include "mrdft/Factory.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

bool initial_guess::nao::setup(OrbitalVector &Phi, double prec, const Nuclei &nucs, int n_mix, double alpha_mix) {
    if (Phi.size() == 0) return false;

    auto restricted = (orbital::size_singly(Phi)) ? false : true;
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Compute initial orbitals");
    print_utils::text(0, "Method      ", "Diagonalize NAO Hamiltonian with optional mixing");
    print_utils::text(0, "Precision   ", print_utils::dbl_to_str(prec, 5, true));
    print_utils::text(0, "Restricted  ", (restricted) ? "True" : "False");
    print_utils::text(0, "Functional  ", "LDA (SVWN5)");
    print_utils::text(0, "AO basis    ", "3-21G");
    mrcpp::print::separator(0, '~', 2);

    bool use_pp = false;
    for (int i = 0; i < nucs.size(); i++) {
        if (nucs[i].hasPseudopotential()) {
            use_pp = true;
            break;
        }
    }

    // Make Fock operator contributions
    Timer t_tot, t_lap;
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);

    mrdft::Factory xc_factory(*MRA);
    xc_factory.setSpin(false);
    xc_factory.setFunctional("SLATERX", 1.0);
    xc_factory.setFunctional("VWN5C", 1.0);
    auto mrdft_p = xc_factory.build();
    MomentumOperator p(D_p);
    CoulombOperator J(P_p);
    XCOperator XC(mrdft_p);
    RankZeroOperator V = J + XC;

    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "NAO Initial Guess");
    if (plevel == 1) mrcpp::print::time(1, "Initializing operators", t_lap);

    // Compute Coulomb density
    t_lap.start();
    Density &rho_j = J.getDensity();

    Density rho_new(false);
    initial_guess::nao::project_atomic_densities(prec, rho_new, nucs);
    double total_charge = 0.0;
    for (int i = 0; i < Phi.size(); i++) {
        total_charge += Phi[i].occ();
    }
    double rho_int = rho_new.integrate().real();
    if (std::abs(rho_int - total_charge) > 1e-6) {
        std::cout << "Total charge: " << total_charge << std::endl;
        std::cout << "Integrated charge: " << rho_int << std::endl;
        rho_new.rescale(total_charge / rho_int);
    }


    // Compute XC density
    Density &rho_xc = XC.getDensity(DensityType::Total);
    // mrcpp::cplxfunc::deep_copy(rho_j, rho_j);

    std::shared_ptr<NuclearOperator> V_nuc;
    std::shared_ptr<ProjectorOperator> P;

    if (use_pp) {
        XC.setNuclei(std::make_shared<Nuclei>(nucs));
        Nuclei nucs_pp;
        Nuclei nucs_all_el;
        for (int i = 0; i < nucs.size(); i++) {
            if (nucs[i].hasPseudopotential()) {
                nucs_pp.push_back(nucs[i]);
            } else {
                nucs_all_el.push_back(nucs[i]);
            }
        }
        P = std::make_shared<ProjectorOperator>(nucs_pp, prec);
        std::string model_pp = "point_like";
        // NuclearOperator V_nuc_all_el(nucs_all_el, prec, prec, false, model_pp);
        V_nuc = std::make_shared<NuclearOperator>(nucs_all_el, prec, prec, false, model_pp);
        model_pp = "pp";
        NuclearOperator pp_nuc(nucs_pp, prec, prec, false,  model_pp);
        V_nuc->add(pp_nuc);
        V = V + (*P);
        // V_nuc_ptr = std::make_shared<NuclearOperator>(V_nuc_all_el);
    } else{
        // NuclearOperator V_nuc(nucs, prec);
        // std::make_shared<NuclearOperator>(V_nuc);
        V_nuc = std::make_shared<NuclearOperator>(nucs, prec);
        // V = V;
    }

    if (plevel == 1) mrcpp::print::time(1, "Projecting GTO density", t_lap);

    t_lap.start();
    OrbitalVector Psi;
    initial_guess::nao::project_atomic_orbitals(prec, Psi, nucs);

    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);


    if (plevel == 1) mrcpp::print::time(1, "Projecting GTO AOs", t_lap);
    if (plevel == 2) mrcpp::print::header(2, "Building Fock operator");
    t_lap.start();
    if (plevel == 2) mrcpp::print::footer(2, t_lap, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_lap);

    p.setup(prec);
    // t_tilde stor
    ComplexMatrix kin_and_nuc_mat = qmoperator::calc_kinetic_matrix(p, Psi, Psi);

    RankZeroOperator V_nuc_op = *V_nuc;
    V_nuc_op.setup(prec);
    kin_and_nuc_mat = kin_and_nuc_mat + V_nuc_op(Psi, Psi);

    for (int imix = 0; imix < n_mix; imix++) {
        std::cout << "Mixing iteration " << imix << std::endl;
        // Compute Fock matrix
        mrcpp::print::header(2, "Diagonalizing Fock matrix");

        mrcpp::cplxfunc::deep_copy(rho_j, rho_new);
        mrcpp::cplxfunc::deep_copy(rho_xc, rho_new);

        V.setup(prec);

        if (imix > 0) {
            double e_kin = qmoperator::calc_kinetic_trace(p, Phi);
            double e_en = V_nuc->trace(Phi).real();
            double e_ee = 0.5 * J.trace(Phi).real();
            double e_xc = XC.getEnergy();
            // double e_pot = V.trace(Phi).real();
            double e_nl = 0.0;
            if (use_pp) {
                e_nl = P->trace(Phi).real();
            }
            double e_nuc = chemistry::compute_nuclear_repulsion(nucs);
            double e_tot = e_kin + e_en + e_ee + e_xc + e_nuc + e_nl;
            std::cout << "Initial LDA energy: " << e_tot << std::endl;
        }

        // ComplexMatrix U = initial_guess::core::diagonalize(Psi, kin_and_nuc_mat, V);
        Timer t1;
        ComplexMatrix f_tilde = kin_and_nuc_mat + V(Psi, Psi);
        ComplexMatrix f = S_m12.adjoint() * f_tilde * S_m12;
        mrcpp::print::separator(2, '-');
        mrcpp::print::time(1, "Computing NAO Fock matrix", t1);

        DoubleVector eig;
        ComplexMatrix U = math_utils::diagonalize_hermitian_matrix(f, eig);
        mrcpp::print::time(1, "Diagonalizing NAO Fock matrix", t1);
        U = S_m12 * U;


        auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
        auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);
        initial_guess::core::rotate_orbitals(Phi, prec, U, Psi);
        initial_guess::core::rotate_orbitals(Phi_a, prec, U, Psi);
        initial_guess::core::rotate_orbitals(Phi_b, prec, U, Psi);
        Phi = orbital::adjoin(Phi, Phi_a);
        Phi = orbital::adjoin(Phi, Phi_b);

        density::compute(prec, rho_new, Phi, DensityType::Total);
        rho_new.rescale(1 - alpha_mix);
        rho_new.add(alpha_mix, rho_j);
        V.clear();
    }

    p.clear();
    V_nuc_op.clear();

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
    return true;
}

void initial_guess::nao::project_atomic_densities(double prec, Density &rho, const Nuclei &nucs){
    std::string source_dir = HIRSHFELD_SOURCE_DIR;
    std::string install_dir = HIRSHFELD_INSTALL_DIR;
    std::string data_dir = "";
    // check if data_dir exists
    if (std::filesystem::exists(install_dir)) {
        data_dir = install_dir + "/lda/";
    } else if (std::filesystem::exists(source_dir)) {
        data_dir = source_dir + "/lda/";
    } else {
        MSG_ABORT("Hirshfeld data directory not found");
    }
    std::vector<interpolation_utils::PolyInterpolator> atomic_densities;

    Density rho_loc(false);
    for (const auto &nuc : nucs) {
        std::string element = nuc.getElement().getSymbol();
        std::string file = data_dir + "/" + element + ".density";

        Eigen::VectorXd rhoGrid, rGrid;
        mrchem::density::readAtomicDensity(file, rGrid, rhoGrid);

        interpolation_utils::PolyInterpolator atomic_density(rGrid, rhoGrid);
        atomic_densities.push_back(atomic_density);
    }

    for (int i = 0; i < nucs.size(); i++) {

        if (mrcpp::my_orb(i)) {
            bool need_rescale = nucs[i].getCharge() != nucs[i].getAtomicNumber();
            mrcpp::ComplexFunction atomic_density_mw;

            if (need_rescale) {
                double cm = nucs[i].getPseudopotentialData()->getMaxRpp() * .5;
                double k = 10;
                double temp = - (1.0 / (1.0 + std::exp(cm * k)));
                mrcpp::Coord<3> nucPos = nucs[i].getCoord();
                auto rho_analytic = [atomic_densities, nucPos, cm, k, i, temp](const mrcpp::Coord<3> &r) {
                    double rr = std::sqrt((r[0] - nucPos[0]) * (r[0] - nucPos[0])
                        + (r[1] - nucPos[1]) * (r[1] - nucPos[1])
                        + (r[2] - nucPos[2]) * (r[2] - nucPos[2]));
                    return atomic_densities[i].evalfLeftNoRightZero(rr) * ((1.0 / (1.0 + std::exp(-k * (rr - cm)))) + temp);
                };
                mrcpp::cplxfunc::project(atomic_density_mw, rho_analytic, mrcpp::NUMBER::Real, prec);
                ComplexDouble integral = atomic_density_mw.integrate();
                double rescale = nucs[i].getCharge() / integral.real();
                // std::cout << "rescale: " << rescale << " integral: " << integral.real() << std::endl;
                atomic_density_mw.rescale(rescale);
            } else {
                mrcpp::Coord<3> nucPos = nucs[i].getCoord();
                auto rho_analytic = [atomic_densities, nucPos, i](const mrcpp::Coord<3> &r) {
                    double rr = std::sqrt((r[0] - nucPos[0]) * (r[0] - nucPos[0])
                        + (r[1] - nucPos[1]) * (r[1] - nucPos[1])
                        + (r[2] - nucPos[2]) * (r[2] - nucPos[2]));
                    return atomic_densities[i].evalfLeftNoRightZero(rr);
                };
                mrcpp::cplxfunc::project(atomic_density_mw, rho_analytic, mrcpp::NUMBER::Real, prec);
            }
            rho_loc.add(1.0, atomic_density_mw);
        }
    }
    density::allreduce_density(prec, rho, rho_loc);
}

class AnalyticOrbital : public mrcpp::RepresentableFunction<3> {
public:
    AnalyticOrbital(int l, int m, const mrcpp::Coord<3> &pos, interpolation_utils::PolyInterpolator rad) {
        this->l = l;
        this->m = m;
        this->pos = pos;
        this->spherical_harmonic = get_spherical_harmonics(l, m);
        this->radial_function = std::make_shared<interpolation_utils::PolyInterpolator>(rad);
    }

    double evalf(const mrcpp::Coord<3> &r) const override {
        double rnorm = std::sqrt((r[0] - pos[0]) * (r[0] - pos[0])
            + (r[1] - pos[1]) * (r[1] - pos[1])
            + (r[2] - pos[2]) * (r[2] - pos[2]));
        mrcpp::Coord<3> rp = {r[0] - pos[0], r[1] - pos[1], r[2] - pos[2]};
        mrcpp::Coord<3> rnormp;
        if (rnorm > 1e-10) {
            rnormp = {rp[0] / rnorm, rp[1] / rnorm, rp[2] / rnorm};
        } else {
            rnormp = {1, 0, 0};
        }
        return radial_function->evalfLeftNoRightZero(rnorm) * spherical_harmonic(rnormp, 1.0);
        
    }
    protected:
    bool isVisibleAtScale(int scale, int nQuadPts) const override {
        // double stdDeviation = 1.0;
        // auto visibleScale = static_cast<int>(std::floor(std::log2(nQuadPts * 5.0 * stdDeviation)));
        // std::cout << "visibleScale: " << visibleScale << " scale: " << scale << std::endl;
        int visibleScale = -5;
        return (scale >= visibleScale);
    }

private:
    int l;
    int m;
    mrcpp::Coord<3> pos;
    std::shared_ptr<interpolation_utils::PolyInterpolator> radial_function;

    double (*spherical_harmonic)(const std::array<double, 3> &r, const double &normr);

};

void initial_guess::nao::project_atomic_orbitals(double prec, OrbitalVector &Phi, const Nuclei &nucs) {
    std::string source_dir = HIRSHFELD_SOURCE_DIR;
    std::string install_dir = HIRSHFELD_INSTALL_DIR;
    OrbitalVector phi_beta;
    std::string data_dir = "";
    // check if data_dir exists
    if (std::filesystem::exists(install_dir)) {
        data_dir = install_dir + "/lda/";
    } else if (std::filesystem::exists(source_dir)) {
        data_dir = source_dir + "/lda/";
    } else {
        MSG_ABORT("Hirshfeld data directory not found");
    }

    for (int iNuc = 0; iNuc < nucs.size(); iNuc++) {
        std::string element = nucs[iNuc].getElement().getSymbol();
        std::string file = data_dir + "/" + element + ".json";
        std::ifstream ifs(file);
        nlohmann::json orbs_json = nlohmann::json::parse(ifs);
        ifs.close();
        // std::cout << "Orbitals for " << element << " are: " << orbs_json << std::endl;
        std::vector<double> r_vec = orbs_json["rgrid"];
        Eigen::VectorXd rGrid(r_vec.size());
        for (int i = 0; i < r_vec.size(); i++) {
            rGrid(i) = r_vec[i];
        }
        for (auto it=orbs_json.begin(); it!=orbs_json.end(); it++) {
            if (it.key() == "rgrid") continue;
            // std::cout << "key: " << it.key() << std::endl;
            char ang_mom_char = it.key()[1];
            std::string ang_mom = std::string(1, ang_mom_char);
            int l;
            if (ang_mom == "s") l = 0;
            else if (ang_mom == "p") l = 1;
            else if (ang_mom == "d") l = 2;
            else if (ang_mom == "f") l = 3;
            else if (ang_mom == "g") l = 4;
            else l = -1;
            std::vector<double> coeffs = it.value();
            Eigen::VectorXd coeffVec(coeffs.size());
            // std::cout << "ang_mom: " << ang_mom << " l: " << l << std::endl;
            // std::cout << "coeffs: " << coeffs.size() << std::endl;
            for (int i = 0; i < coeffs.size(); i++) {
                coeffVec(i) = coeffs[i];
            }
            interpolation_utils::PolyInterpolator rad_func(rGrid, coeffVec);
            for (int m = -l; m <= l; m++) {
                mrcpp::Coord<3> pos = nucs[iNuc].getCoord();
                AnalyticOrbital orb(l, m, pos, rad_func);
                Orbital orb_mw;
                mrcpp::cplxfunc::project(orb_mw, orb, mrcpp::NUMBER::Real, prec);
                double nrm = orb_mw.norm();
                orb_mw.rescale(1 / nrm);
                // std::cout << "n nodes " << orb_mw.getNNodes(mrcpp::NUMBER::Total) << std::endl;
                // std::cout << "norm " << orb_mw.norm() << std::endl;
                Phi.push_back(orb_mw);
            }
        } 
    }
}

} // namespace mrchem
