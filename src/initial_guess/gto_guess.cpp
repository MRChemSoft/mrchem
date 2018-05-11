#include "MRCPP/MWFunctions"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "mathutils.h"

#include "gto_guess.h"
#include "gto_utils/OrbitalExp.h"
#include "gto_utils/Intgrl.h"

#include "Molecule.h"
#include "Orbital.h"

using mrcpp::GaussExp;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace gto_guess {

// Forward declare helper functions
void project(double prec, OrbitalVector &Phi, gto_utils::OrbitalExp &gto_exp);

} //namespace gto_guess


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
OrbitalVector gto_guess::initial_guess(double prec,
                                       const Molecule &mol,
                                       const std::string &bas_file,
                                       const std::string &mo_file) {
    Printer::printHeader(0, "Setting up occupied orbitals (closed-shell)");
    println(0, "    n  Spin  Occ                           SquareNorm");
    Printer::printSeparator(0, '-');
    Timer timer;

    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity();   //multiplicity
    int Ne = mol.getNElectrons();       //total electrons
    int Nd = Ne - (mult - 1);           //doubly occupied electrons
    if (Nd%2 != 0) MSG_FATAL("Invalid multiplicity");

    // Read basis set file
    gto_utils::Intgrl intgrl(bas_file);

    // Setup AO basis
    gto_utils::OrbitalExp gto_exp(intgrl);

    // Read MO file and rotate into MO basis
    DoubleMatrix MO = mathutils::read_matrix_file(mo_file);
    gto_exp.rotate(MO.transpose());

    // Setup empty orbitals
    OrbitalVector Phi;
    for (int i = 0; i < Nd/2; i++)  Phi.push_back(SPIN::Paired);

    // Project GTO expansion
    gto_guess::project(prec, Phi, gto_exp);

    timer.stop();
    Printer::printFooter(0, timer, 2);

    return Phi;
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
OrbitalVector gto_guess::initial_guess(double prec,
                                       const Molecule &mol,
                                       const std::string &bas_file,
                                       const std::string &moa_file,
                                       const std::string &mob_file) {
    Printer::printHeader(0, "Setting up occupied orbitals (open-shell)");
    println(0, "    n  Spin  Occ                           SquareNorm");
    Printer::printSeparator(0, '-');
    Timer timer;

    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity();   //multiplicity
    int Ne = mol.getNElectrons();       //total electrons
    int Nd = Ne - (mult - 1);           //paired electrons
    if (Nd%2 != 0) MSG_FATAL("Invalid multiplicity");
    int Na = Nd/2 + (mult - 1);         //alpha electrons
    int Nb = Nd/2;                      //beta electrons

    // Read basis set file
    gto_utils::Intgrl intgrl(bas_file);

    // Alpha orbitals
    OrbitalVector Phi_a;
    {
        // Setup AO basis
        gto_utils::OrbitalExp gto_exp(intgrl);

        // Read MO file and rotate into MO basis
        DoubleMatrix MO_a = mathutils::read_matrix_file(moa_file);
        gto_exp.rotate(MO_a.transpose());

        // Setup empty orbitals
        for (int i = 0; i < Na; i++)  Phi_a.push_back(SPIN::Alpha);

        // Project GTO expansion
        gto_guess::project(prec, Phi_a, gto_exp);
    }

    // Beta orbitals
    OrbitalVector Phi_b;
    {
        // Read MO file
        DoubleMatrix MO_b = mathutils::read_matrix_file(mob_file);

        // Setup AO basis and rotate into MOs
        gto_utils::OrbitalExp gto_exp(intgrl);
        gto_exp.rotate(MO_b.transpose());

        // Setup empty orbitals
        for (int i = 0; i < Nb; i++)  Phi_b.push_back(SPIN::Beta);

        // Project GTO expansion
        gto_guess::project(prec, Phi_b, gto_exp);
    }

    timer.stop();
    Printer::printFooter(0, timer, 2);

    // Collect orbitals into one vector
    return orbital::adjoin(Phi_a, Phi_b);
}

/** @brief Project the N first GTO expansions of the MO basis
 *
 * @param prec: precision used in projection
 * @param Phi: vector or MW orbitals
 * @param gto_exp: vector of GTO orbitals
 *
 * Projects the GTO orbital expansions into the corresponding MW orbitals.
 * The length is decided by the MW vector and it is assumed that the GTO
 * vector is at least of the same size.
 *
 */
void gto_guess::project(double prec, OrbitalVector &Phi, gto_utils::OrbitalExp &gto_exp) {
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            Phi[i].alloc(NUMBER::Real);
            mrcpp::project(prec, Phi[i].real(), gto_exp[i]);
            printout(0, std::setw(5)  << i);
            printout(0, std::setw(5)  << Phi[i].printSpin());
            printout(0, std::setw(5)  << Phi[i].occ());
            printout(0, std::setw(44) << Phi[i].norm() << std::endl);
        }
    }
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

