#ifndef ORBITALVECTOR_H
#define ORBITALVECTOR_H

#pragma GCC system_header
#include <Eigen/Core>
#pragma GCC system_header
#include <Eigen/Eigenvalues>

#include <vector>

#include "Orbital.h"

class OrbitalVector {
public:
    OrbitalVector(int n_orbs);
    OrbitalVector(int n_alpha, int n_beta);
    OrbitalVector(int ne, int mult, bool rest);
    OrbitalVector(const OrbitalVector &orb_set);
    OrbitalVector& operator=(const OrbitalVector &orb_set) { NOT_IMPLEMENTED_ABORT; }
    virtual ~OrbitalVector();

    void deepCopy(OrbitalVector &inp);
    void shallowCopy(const OrbitalVector &inp);

    void push_back(int n_orbs, int occ, int spin);
    void push_back(Orbital& Orb);
    void pop_back(bool free = true);
    void clear(bool free = true);
    void clearVec(bool free = true);

    void normalize();

    Eigen::MatrixXcd calcOverlapMatrix();
    Eigen::MatrixXcd calcOverlapMatrix(OrbitalVector &ket);
    Eigen::MatrixXcd calcOverlapMatrix_P(OrbitalVector &ket);
    Eigen::MatrixXcd calcOverlapMatrix_P_H(OrbitalVector &ket);

    int size() const { return this->orbitals.size(); }
    int getNOccupied() const;
    int getNEmpty() const;
    int getNSingly() const;
    int getNDoubly() const;
    int getNPaired() const;
    int getNAlpha() const;
    int getNBeta() const;
    int getNElectrons(int spin = Orbital::Paired) const;
    int getMultiplicity() const;

    void setSpins(const Eigen::VectorXi &spins);
    void setErrors(const Eigen::VectorXd &errors);
    void setOccupancies(const Eigen::VectorXi &occ);

    Eigen::VectorXi getSpins() const;
    Eigen::VectorXd getNorms() const;
    Eigen::VectorXd getErrors() const;
    Eigen::VectorXd getSquareNorms() const;
    Eigen::VectorXi getOccupancies() const;

    bool isConverged(double prec) const;
    double calcTotalError() const;

    const Orbital *operator[](int i) const { return this->orbitals[i]; }
    Orbital *operator[](int i) { return this->orbitals[i]; }

    const Orbital &getOrbital(int i) const;
    Orbital &getOrbital(int i);

    void replaceOrbital(int i, Orbital **orb);

    int printTreeSizes() const;

    void send_OrbVec(int dest, int tag, int* OrbsIx, int start, int maxcount);
    void Isend_OrbVec(int dest, int tag, int* OrbsIx, int start, int maxcount);
    void Rcv_OrbVec(int source, int tag, int* OrbsIx, int& workOrbVecIx);
    void getOrbVecChunk(int* myOrbsIx, OrbitalVector &rcvOrbs, int* rcvOrbsIx, int size, int& iter0);
    void getOrbVecChunk_sym(int* myOrbsIx, OrbitalVector &rcvOrbs, int* rcvOrbsIx, int size, int& iter0);

    friend std::ostream& operator<<(std::ostream &o, OrbitalVector &orb_set) {
        int oldPrec = TelePrompter::setPrecision(15);
        o << "*OrbitalVector: ";
        o << std::setw(4) << orb_set.size()          << " orbitals  ";
        o << std::setw(4) << orb_set.getNOccupied()  << " occupied  ";
        o << std::setw(4) << orb_set.getNElectrons() << " electrons " << std::endl;
        o << "------------------------------";
        o << "------------------------------\n";
        o << "   n    sqNorm               Occ Spin  Error\n";
        o << "------------------------------";
        o << "------------------------------\n";
        for (int i = 0; i < orb_set.size(); i++) {
            Orbital *orb = orb_set[i];
            o << std::setw(4) << i;
            if (orb != 0) {
                o << *orb;
            } else {
                o << std::endl;
            }
        }
        o << "------------------------------";
        o << "------------------------------\n";
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    //Data
    //LUCA: ... or maybe just have two sets here?
    std::vector<Orbital *> orbitals;
};

#endif // ORBITALVECTOR_H
