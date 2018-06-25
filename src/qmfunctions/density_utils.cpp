#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/math_utils.h"
#include "utils/RRMaximizer.h"

#include "qmfunctions.h"
#include "Density.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace density {
/****************************************
 * Density related standalone functions *
 ****************************************/

Density compute(double prec, Orbital phi, int spin) {
    occ = calc_density_occupancy(phi.spin(), spin, phi.occ());
    Density rho(spin);
    if (std::abs(occ) < mrcpp::MachineZero) {
        rho.real().setZero();
        rho.imag().setZero();
        return 0;
    }

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*real_2, rho.real());
        mrcpp::multiply(prec, *real_2, occ, phi.real(), phi.real());
        sum_vec.push_back(std::make_tuple(1.0, real_2));
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*imag_2, rho.real());
        mrcpp::multiply(prec, *imag_2, occ, phi.imag(), phi.imag());
        sum_vec.push_back(std::make_tuple(-1.0, imag_2));  //@stig: this should be -1.0. It was 1.0...
    }
    mrcpp::build_grid(rho.real(), sum_vec);
    mrcpp::add(-1.0, rho.real(), sum_vec, 0);
    mrcpp::clear(sum_vec, true);
    return rho;
}

/*
Density compute_singlet_tranition_density(double prec, Density &rho, Orbital phi_0, Orbital phi_1, int spin) {

    double occ = calc_density_occupancy(phi_0.spin(), spin, phi_0.occ());

    if (std::abs(occ) < mrcpp::MachineZero) {
        rho.real().setZero();
        rho.imag().setZero();
        return;
    }

    qmfunction::multiply(phi_0, phi_0.conjugate(), phi_1, phi_1.conjugate(), dens, occ, prec);

}
*/

Density compute(double prec, Density &rho, OrbitalVector &Phi, int spin) {
    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i

    DensityVector dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        Density *rho_i = new Density(*MRA);
        mrcpp::copy_grid(*rho_i, rho);
        density::compute(mult_prec, *rho_i, Phi[i], spin);
        dens_vec.push_back(std::make_tuple(1.0, rho_i));
    }

    // Adaptive prec addition if more than 5 contributions,
    // otherwise addition on union grid
    if (dens_vec.size() > 5 and add_prec > 0.0) {  
        mrcpp::add(add_prec, rho, dens_vec);
    } else if (dens_vec.size() > 0) {
        mrcpp::build_grid(rho, dens_vec); //LUCA: this does not seem to fit with XCPotential (grid kept as provided by MRDFT)
        mrcpp::add(-1.0, rho, dens_vec, 0);
    }
    mrcpp::clear(dens_vec, true);
}

/*
Density compute(double prec, Density &rho, OrbitalVector &Phi_0,
                      OrbitalVector &Phi_1, int spin) {
    double mult_prec = prec;            // prec for \rho_i = \phi_0_i * \phi_1_i
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i
    
    DensityVector dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        Density *rho_i = new Density();
        mrcpp::copy_grid(*rho_i, rho);
        density::compute(mult_prec, *rho_i, Phi_0[i], Phi_1[i], spin);
        dens_vec.push_back(std::make_tuple(1.0, rho_i));
    }

    // Adaptive prec addition if more than 5 contributions,
    // otherwise addition on union grid
    if (dens_vec.size() > 5 and add_prec > 0.0) {
        mrcpp::add(add_prec, rho, dens_vec);
    } else if (dens_vec.size() > 0) {
        mrcpp::build_grid(rho, dens_vec);
        mrcpp::add(-1.0, rho, dens_vec, 0);
    }
    mrcpp::clear(dens_vec, true);

    bool ac = (phi_0.hasReal() and phi_1.hasReal());
    bool bd = (phi_0.hasImag() and phi_1.hasImag());
    bool ad = (phi_0.hasReal() and phi_1.hasImag());
    bool bc = (phi_0.hasImag() and phi_1.hasReal());
    double conj_0 = (phi_0.conjugate()) ? -1.0 : 1.0;
    double conj_1 = (phi_1.conjugate()) ? -1.0 : 1.0;

    if (ac)                   multiply(phi_0.real(), phi_1.real());
    if (bd) conj_0 * conj_1 * multiply(phi_0.imag(), phi_1.imag());
    if (ad) conj_1 *          multiply(phi_0.real(), phi_1.imag());
    if (bc) conj_0 *          multiply(phi_0.imag(), phi_1.real());
    }
}
*/

/** @brief computes the occupancy of a given density term
 *
 * parma[in] spin_orb the spin of the corresponding orbital 
 * parma[in] spin_dens the spin of the density term
 * parma[in] occ_orb the occupancy value of the oribital
 *
 * Note: this is valid for a density of a single Slater determinant:
 * both close and open shell
 */
double density::calc_density_occupancy(int spin_orb, int spin_dens, double occ_orb) {

    Double occ_a(0.0), occ_b(0.0), occ_p(0.0);
    switch (spin_orb){
      case (SPIN::Alpha):  occ_a = (double) occ_orb; break;
      case (SPIN::Beta):   occ_b = (double) occ_orb; break;
      case (SPIN::Paired): occ_p = (double) occ_orb; break;
      case default: MSG_ABORT("Invalid orbital spin");
    }
    
    double occ(0.0);
    switch (spin_dens) {
      case (DENSITY::Total): occ = occ_a + occ_b + occ_p; break;
      case (DENSITY::Alpha): occ = occ_a + 0.5*occ_p;     break;
      case (DENSITY::Beta):  occ = occ_b + 0.5*occ_p;     break;
      case (DENSITY::Spin):  occ = occ_a - occ_b;         break;
      case default: MSG_ABORT("Invalid density spin");   
    }
    return occ;
}

/** @brief computes the occupancy of a given transition density term
 *
 * parma[in] spin_orb the spin of the corresponding orbital 
 * parma[in] spin_dens the spin of the density term
 * parma[in] occ_orb the occupancy value of the oribital
 *
 * Note: this is valid for a density of a single Slater determinant:
 * both close and open shell
 */
double density::calc_singlet_density_occupancy(spin_orb_0, spin_orb_1, spin_dens, occ_orb_0, occ_orb_1) {

    if(spin_orb_0 != spin_orb_1) MSG_ABORT("Inconsistent spins"); 
    if(occ_orb_0  != occ_orb_1)  MSG_ABORT("Inconsistent occupancies"); 
    return calc_density_occupancy(spin_orb_0, spin_dens, occ_orb_0;
}

        
double density::calc_triplet_density_occupancy(spin_orb_0, spin_orb_1, spin_dens, occ_orb_0, occ_orb_1) {
        NOT_IMPLEMENTED_ABORT;
    }

}

} //namespace density

} //namespace mrchem
