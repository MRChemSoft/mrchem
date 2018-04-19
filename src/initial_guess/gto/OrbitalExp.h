#pragma once

#include <vector>

#include "MRCPP/Gaussians"

#include "mrchem.h"

namespace mrchem {

namespace gto_guess {
class Intgrl;

class OrbitalExp final {
public:
    OrbitalExp(Intgrl &intgrl);
    ~OrbitalExp();

    int size() const { return this->orbitals.size(); }
    int getAngularMomentum(int n) const;

    mrcpp::GaussExp<3> &operator[](int n) { return *this->orbitals[n]; }

    void rotate(const DoubleMatrix &U);

protected:
    bool cartesian;
    std::vector<mrcpp::GaussExp<3> *> orbitals;

    void readAOExpansion(Intgrl &intgrl);
    void transformToSpherical();
};

} //namespace gto_guess
} //namespace mrchem
