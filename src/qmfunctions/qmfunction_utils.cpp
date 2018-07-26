#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/math_utils.h"
#include "utils/RRMaximizer.h"

#include "qmfunction_utils.h"
#include "Orbital.h"
#include "OrbitalIterator.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA


/****************************************
 * Density related standalone functions *
 ****************************************/

void density::compute(double prec, Density &rho, Orbital phi, int spin) {
    double occ_a(0.0), occ_b(0.0), occ_p(0.0);
    if (phi.spin() == SPIN::Alpha)  occ_a = (double) phi.occ();
    if (phi.spin() == SPIN::Beta)   occ_b = (double) phi.occ();
    if (phi.spin() == SPIN::Paired) occ_p = (double) phi.occ();

    double occ(0.0);
    if (spin == DENSITY::Total) occ = occ_a + occ_b + occ_p;
    if (spin == DENSITY::Alpha) occ = occ_a + 0.5*occ_p;
    if (spin == DENSITY::Beta)  occ = occ_b + 0.5*occ_p;
    if (spin == DENSITY::Spin)  occ = occ_a - occ_b;

    if (std::abs(occ) < mrcpp::MachineZero) {
        rho.setZero();
        return;
    }

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*real_2, phi.real());
        mrcpp::square(prec, *real_2, phi.real());
        sum_vec.push_back(std::make_tuple(occ, real_2));
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*imag_2, phi.imag());
        mrcpp::square(prec, *imag_2, phi.imag());
        sum_vec.push_back(std::make_tuple(occ, imag_2));
    }
    mrcpp::build_grid(rho, sum_vec);
    mrcpp::add(-1.0, rho, sum_vec, 0);
    mrcpp::clear(sum_vec, true);
}

void density::compute(double prec, Density &rho, OrbitalVector &Phi, int spin) {
    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i

    DensityVector dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            Density *rho_i = new Density(*MRA);
            mrcpp::copy_grid(*rho_i, rho);
            density::compute(mult_prec, *rho_i, Phi[i], spin);
            dens_vec.push_back(std::make_tuple(1.0, rho_i));
        }
    }

    if (add_prec > 0.0) {
        mrcpp::add(add_prec, rho, dens_vec);
    } else if (dens_vec.size() > 0) {
        mrcpp::build_grid(rho, dens_vec);
        mrcpp::add(-1.0, rho, dens_vec, 0);
    }
    mrcpp::clear(dens_vec, true);

    mpi::reduce_density(rho, mpi::comm_orb);
    mpi::broadcast_density(rho, mpi::comm_orb);
}

} //namespace mrchem
