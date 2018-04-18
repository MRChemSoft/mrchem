#include <Eigen/Eigenvalues>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "hydrogen_guess.h"

#include "HydrogenFunction.h"
#include "Molecule.h"
#include "Nucleus.h"
#include "Orbital.h"

#include "NuclearOperator.h"
#include "KineticOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** Helper struct to get the orbital ordering right
 *
 *  First index energy level (n)
 *  Second index angular momentum (l)
 */
namespace hydrogen_guess {
int PT[29][2] = {
   /*s*/
   {1,0},                  /*p*/
   {2,0},                  {2,1},
   {3,0},            /*d*/ {3,1},
   {4,0},            {3,2},{4,1},
   {5,0},      /*f*/ {4,2},{5,1},
   {6,0},      {4,3},{5,2},{6,1},
   {7,0},/*g*/ {5,3},{6,2},{7,1},
   {8,0},{5,4},{6,3},{7,2},{8,1},
   {9,0},{6,4},{7,3},{8,2},{9,1}
};

OrbitalVector initial_guess(double prec,
                            const Molecule &mol,
                            bool restricted,
                            int zeta) {
    int mult = mol.getMultiplicity();   //multiplicity
    int Ne = mol.getNElectrons();       //total electrons
    int Nd = Ne - (mult - 1);           //doubly occupied
    if (Nd%2 != 0) MSG_FATAL("Invalid multiplicity");
    if (not restricted) NOT_IMPLEMENTED_ABORT;

    //project AO basis of hydrogen functions
    OrbitalVector Phi = hydrogen_guess::project(prec, mol.getNuclei(), zeta);

    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);

    // Compute core Hamiltonian matrix
    mrcpp::ABGVOperator<3> D(*MRA, 0.5, 0.5);
    KineticOperator T(D);
    NuclearOperator V_nuc(mol.getNuclei(), prec);
    T.setup(prec);
    V_nuc.setup(prec);
    ComplexMatrix t = T(Phi, Phi);
    ComplexMatrix v = V_nuc(Phi, Phi);
    V_nuc.clear();
    T.clear();

    ComplexMatrix F = S_m12.transpose()*(t + v)*S_m12;

    // Diagonalize Hamiltonian matrix
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(F.cols());
    es.compute(F);
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexMatrix U = ei_vec.transpose() * S_m12;

    // Rotate orbitals and fill electrons by Aufbau
    OrbitalVector Psi;
    for (int i = 0; i < Nd/2; i++) {
        ComplexVector v_i = U.row(i);
        Orbital psi_i = orbital::multiply(v_i, Phi, prec);
        Psi.push_back(psi_i);
    }
    orbital::free(Phi);

    if (restricted) {
        if (mult != 1) MSG_FATAL("Restricted open-shell not available");

        //set spin and occupation number
        hydrogen_guess::populate(Psi, Nd, SPIN::Paired);
    } else {
        NOT_IMPLEMENTED_ABORT;
        /*
        OrbitalVector Phi_a = Phi;
        OrbitalVector Phi_b = orbital::deep_copy(Phi);

        int Na = Nd/2 + (mult - 1);     //alpha electrons
        int Nb = Nd/2;                  //beta electrons

        //set spin and occupation number
        hydrogen_guess::populate(Phi_a, Na, SPIN::Alpha);
        hydrogen_guess::populate(Phi_b, Nb, SPIN::Beta);

        Phi.clear();
        Phi = orbital::adjoin(Phi_a, Phi_b);
        */
    }
    return Psi;
}

OrbitalVector project(double prec, const Nuclei &nucs, int zeta) {
    Printer::printHeader(0, "Setting up occupied orbitals");
    println(0, "    N    Atom   Label                     SquareNorm");
    Printer::printSeparator(0, '-');

    Timer timer;
    OrbitalVector Phi;

    const char label[10] = "spdfg";

    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        int minAO = std::ceil(nuc.getElement().getZ()/2.0);
        double Z = nuc.getCharge();
        const double *R = nuc.getCoord();

        int nAO = 0;
        int nShell = 0;
        int zetaReached = 0;
        bool minAOReached = false;
        while (true) {
            int n = hydrogen_guess::PT[nShell][0];
            int l = hydrogen_guess::PT[nShell][1];
            int M = 2*l + 1;

            if (minAOReached and l == 0) zetaReached++;
            if (zetaReached >= zeta) break;

            for (int m = 0; m < M; m++) {
                HydrogenFunction h_func(n, l, m, Z, R);

                Orbital phi_i;
                phi_i.alloc(NUMBER::Real);
                mrcpp::project(prec, phi_i.real(), h_func);

                printout(0, std::setw(5)  << Phi.size());
                printout(0, std::setw(6)  << nuc.getElement().getSymbol() << i+1);
                printout(0, std::setw(6)  << n << label[l]);
                printout(0, std::setw(40) << phi_i.squaredNorm());
                printout(0, std::endl);

                Phi.push_back(phi_i);
                if (++nAO >= minAO) minAOReached = true;
            }
            nShell++;
        }
    }
    timer.stop();
    Printer::printFooter(0, timer, 2);
    return Phi;
}

void populate(OrbitalVector &vec, int N, int spin) {
    int occ = 0;
    if (spin == SPIN::Paired) occ = 2;
    if (spin == SPIN::Alpha) occ = 1;
    if (spin == SPIN::Beta) occ = 1;
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setSpin(spin);
        if (i < N) {
            vec[i].setOcc(occ);
        } else {
            vec[i].setOcc(0);
        }
    }
}

} //namespace hydrogen_guess

} //namespace mrchem
