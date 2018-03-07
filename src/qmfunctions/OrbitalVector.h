#pragma once

#include <vector>
using std::vector;

#include "Orbital.h"

namespace mrchem {

class OrbitalVector {
public:
    OrbitalVector() : in_use(false) { }
    OrbitalVector(int ne, int mult, bool rest);
    ~OrbitalVector() { }

    OrbitalVector(const OrbitalVector &inp_vec);
    OrbitalVector& operator=(const OrbitalVector &inp_vec);

    OrbitalVector deepCopy();
    OrbitalVector paramCopy() const;

    void clear();
    void free();

    void push_back(Orbital inp) { this->orbitals.push_back(inp); }

    void adjoin(OrbitalVector &inp_b);
    OrbitalVector disjoin(int spin);

    void normalize();
    void orthogonalize();
    void orthogonalize(OrbitalVector inp_vec);

    void rotate(const ComplexMatrix &U, double prec = -1.0);

    int size() const { return this->orbitals.size(); }
    int getNOccupied() const;
    int getNEmpty() const;
    int getNSingly() const;
    int getNDoubly() const;
    int getNPaired() const;
    int getNAlpha() const;
    int getNBeta() const;
    int getNElectrons(int spin = SPIN::Paired) const;
    int getMultiplicity() const;

    void setSpins(const IntVector &spins);
    void setErrors(const DoubleVector &errors);
    void setOccupancies(const IntVector &occ);

    IntVector getSpins() const;
    IntVector getOccupancies() const;
    DoubleVector getNorms() const;
    DoubleVector getErrors() const;
    DoubleVector getSquaredNorms() const;

    const Orbital& operator[](int i) const { return this->orbitals[i]; }
    Orbital& operator[](int i) { return this->orbitals[i]; }

    friend std::ostream& operator<<(std::ostream &o, OrbitalVector orb_vec) { return orb_vec.print(o); }

protected:
    bool in_use;
    std::vector<Orbital> orbitals;

    std::ostream& print(std::ostream &o) const;
};

} //namespace mrchem
