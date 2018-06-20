#include "MRCPP/MWFunctions"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "math_utils.h"

#include "initial_guess/gto.h"
#include "gto_utils/OrbitalExp.h"
#include "gto_utils/Intgrl.h"

#include "Molecule.h"
#include "Orbital.h"

using mrcpp::GaussExp;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Produce an initial guess of orbitals
 *
 * @param prec: precision used in projection
 * @param mol: molecule
 * @param bas_file: basis set file (LSDalton format)
 * @param mo_file: file with MO coefficients
 *
 * Sets up a precomputed MO basis from a spin restricted LSDalton calculation.
 * Requires the LSDalton basis file and the corresponding MO matrix (not in any
 * official format!). The MO file should start with one entry giving the number
 * of AOs, followed by the columns of the MO matrix concatenated into a single
 * column.
 *
 * Projects only the occupied orbitals.
 *
 */
OrbitalVector initial_guess::gto::setup(double prec,
                                        const Molecule &mol,
                                        const std::string &bas_file,
                                        const std::string &mo_file) {
    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity();   //multiplicity
    int Ne = mol.getNElectrons();       //total electrons
    int Nd = Ne - (mult - 1);           //doubly occupied electrons
    if (Nd%2 != 0) MSG_FATAL("Invalid multiplicity");
    int Np = Nd/2;                      //paired orbitals

    // Project GTO expansion
    return initial_guess::gto::project_mo(prec, bas_file, mo_file, SPIN::Paired, 2, Np);
}

/** @brief Produce an initial guess of orbitals
 *
 * @param prec: precision used in projection
 * @param mol: molecule
 * @param bas_file: basis set file (LSDalton format)
 * @param moa_file: file with alpha MO coefficients
 * @param mob_file: file with beta MO coefficients
 *
 * Sets up a precomputed MO basis from an unrestricted LSDalton calculation.
 * Requires the LSDalton basis file and the corresponding MO matrices (not in
 * any official format!). The MO files should start with one entry giving the
 * number of AOs, followed by the columns of the MO matrix concatenated into
 * a single column.
 *
 * Projects only the occupied orbitals of each spin.
 *
 */
OrbitalVector initial_guess::gto::setup(double prec,
                                        const Molecule &mol,
                                        const std::string &bas_file,
                                        const std::string &moa_file,
                                        const std::string &mob_file) {
    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity();   //multiplicity
    int Ne = mol.getNElectrons();       //total electrons
    int Nd = Ne - (mult - 1);           //paired electrons
    if (Nd%2 != 0) MSG_FATAL("Invalid multiplicity");
    int Na = Nd/2 + (mult - 1);         //alpha orbitals
    int Nb = Nd/2;                      //beta orbitals

    // Project orbitals
    OrbitalVector Phi_a = initial_guess::gto::project_mo(prec, bas_file, moa_file, SPIN::Alpha, 1, Na);
    OrbitalVector Phi_b = initial_guess::gto::project_mo(prec, bas_file, mob_file, SPIN::Beta, 1, Nb);

    // Collect orbitals into one vector
    return orbital::adjoin(Phi_a, Phi_b);
}

/** @brief Project the N first GTO expansions of the MO basis
 *
 * @param prec Precision used in projection
 * @param bas_file String containing basis set file
 * @param mo_file String containing MO matrix file
 * @param N Number of orbitals to project
 * @param spin Spin parameter of orbitals
 *
 * @returns Phi: vector or MW orbitals
 *
 * Projects the N first rows of the MO matrix from GTO orbitals into
 * corresponding MW orbitals. All orbitals get the same spin parameter.
 *
 */
OrbitalVector initial_guess::gto::project_mo(double prec,
                                             const std::string &bas_file,
                                             const std::string &mo_file,
                                             int spin,
                                             int occ,
                                             int N) {
    Printer::printHeader(0, "Setting up Gaussian-type MOs");
    println(0, "    n  Spin  Occ                           SquareNorm");
    Printer::printSeparator(0, '-');
    Timer timer;

    // Setup AO basis
    gto_utils::Intgrl intgrl(bas_file);
    gto_utils::OrbitalExp gto_exp(intgrl);
    if (N < 0) N = gto_exp.size();

    // Read MO file (transpose)
    DoubleMatrix MO = math_utils::read_matrix_file(mo_file);
    if (MO.cols() < N) MSG_FATAL("Size mismatch");

    OrbitalVector Phi;
    for (int i = 0; i < N; i++) {
        Orbital phi_i(spin, occ, i%mpi::orb_size);
        if (mpi::my_orb(phi_i)) {
            phi_i.alloc(NUMBER::Real);
            GaussExp<3> mo_i = gto_exp.getMO(i, MO.transpose());
            mrcpp::project(prec, phi_i.real(), mo_i);
        }
        printout(0, std::setw(5)  << i);
        printout(0, std::setw(5)  << phi_i.printSpin());
        printout(0, std::setw(5)  << phi_i.occ());
        printout(0, std::setw(44) << phi_i.norm() << std::endl);
        Phi.push_back(phi_i);
    }

    MPI_Barrier(mpi::comm_orb);
    timer.stop();
    Printer::printFooter(0, timer, 2);

    return Phi;
}

/** @brief Project the N first GTO expansions of the AO basis
 *
 * @param prec Precision used in projection
 * @param bas_file String containing basis set file
 * @param N Number of orbitals to project
 * @param spin Spin parameter of orbitals
 *
 * @returns Phi: vector or MW orbitals
 *
 * Projects the N first Gaussian-type AOs into corresponding MW orbitals.
 * All orbitals get the same spin parameter.
 *
 */
OrbitalVector initial_guess::gto::project_ao(double prec,
                                             const std::string &bas_file,
                                             int spin,
                                             int occ,
                                             int N) {
    Printer::printHeader(0, "Setting up Gaussian-type AOs");
    println(0, "    n  Spin  Occ                           SquareNorm");
    Printer::printSeparator(0, '-');
    Timer timer;

    // Setup AO basis
    gto_utils::Intgrl intgrl(bas_file);
    gto_utils::OrbitalExp gto_exp(intgrl);
    if (N < 0) N = gto_exp.size();

    OrbitalVector Phi;
    for (int i = 0; i < N; i++) {
        Orbital phi_i(spin, occ);
        phi_i.alloc(NUMBER::Real);
        GaussExp<3> ao_i = gto_exp.getAO(i);
        mrcpp::project(prec, phi_i.real(), ao_i);
        printout(0, std::setw(5)  << i);
        printout(0, std::setw(5)  << phi_i.printSpin());
        printout(0, std::setw(5)  << phi_i.occ());
        printout(0, std::setw(44) << phi_i.norm() << std::endl);
        Phi.push_back(phi_i);
    }

    MPI_Barrier(mpi::comm_orb);
    timer.stop();
    Printer::printFooter(0, timer, 2);

    return Phi;
}

/** @brief Project the N first GTO expansions of the MO basis
 *
 * @param prec Precision used in projection
 * @param bas_file String containing basis set file
 * @param mo_file String containing MO matrix file
 * @param N Number of orbitals to project
 * @param spin Spin parameter of orbitals
 *
 * @returns Phi: vector or MW orbitals
 *
 * Projects the N first rows of the MO matrix from GTO orbitals into
 * corresponding MW orbitals. All orbitals get the same spin parameter.
 *
 */
Density* initial_guess::gto::project_density(double prec,
                                             const Nucleus &nuc,
                                             const std::string &bas_file,
                                             const std::string &dens_file) {
    // Setup AO basis
    gto_utils::Intgrl intgrl(bas_file);
    intgrl.getNucleus(0).setCoord(nuc.getCoord());
    gto_utils::OrbitalExp gto_exp(intgrl);

    // Read density matrix file
    DoubleMatrix D = math_utils::read_matrix_file(dens_file);
    GaussExp<3> dens_exp = gto_exp.getDens(D);

    Density *rho = new Density(*MRA);
    mrcpp::project(prec, *rho, dens_exp);
    return rho;
}

} //namespace mrchem


//void OrbitalVector::readVirtuals(const string &bf, const string &mo, int n_occ) {
//    Timer timer;
//    int oldPrec = Printer::setPrecision(15);
//    printout(0, "\n\n=============== Setting up virtual orbitals ");
//    printout(0, "================\n\n");

//    OrbitalExp *moExp = readOrbitalExpansion(bf, mo);
//    for (int a = n_occ; a < moExp->size(); a++) {
//    GaussExp<3> &gtOrb = moExp->getOrbital(a);
//        Orbital *orb_a = new Orbital(2, Orbital::Paired);
//    orb_a->projectFunction(gtOrb);
//        printout(0, "Orbital " << setw(3) << a);
//        println(0, " squareNorm: " << setw(36) << orb_a->getSquareNorm());
//        this->orbitals.push_back(orb_a);
//    }
//    delete moExp;
//    Printer::setPrecision(5);
//    printout(0, "\n================ Elapsed time: ");
//    println(0, timer.elapsed() << " =================\n");
//    Printer::setPrecision(oldPrec);
//}

