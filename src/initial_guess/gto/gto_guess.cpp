#include "MRCPP/MWFunctions"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "gto_guess.h"
#include "OrbitalExp.h"
#include "Intgrl.h"

#include "utils/mathutils.h"
#include "parallel.h"
#include "Molecule.h"
#include "Orbital.h"

using mrcpp::GaussExp;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

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
    Intgrl intgrl(bas_file);

    // Setup AO basis
    gto_guess::OrbitalExp gto_exp(intgrl);

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
    Intgrl intgrl(bas_file);

    // Alpha orbitals
    OrbitalVector Phi_a;
    {
        // Setup AO basis
        OrbitalExp gto_exp(intgrl);

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
        gto_guess::OrbitalExp gto_exp(intgrl);
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


void gto_guess::project(double prec, OrbitalVector &Phi, gto_guess::OrbitalExp &gto_exp) {
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

