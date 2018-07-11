#pragma once

#include <vector>

#include "MRCPP/Gaussians"

#include "mrchem.h"

namespace mrchem {

namespace gto_utils {
class Intgrl;

class OrbitalExp final {
public:
    OrbitalExp(Intgrl &intgrl);
    ~OrbitalExp();

    int size() const { return this->orbitals.size(); }
    int getAngularMomentum(int n) const;

    mrcpp::GaussExp<3> getAO(int i) const { return *this->orbitals[i]; }
    mrcpp::GaussExp<3> getMO(int i, const DoubleMatrix &M) const;
    mrcpp::GaussExp<3> getDens(const DoubleMatrix &D) const;

    void rotate(const DoubleMatrix &U);

protected:
    bool cartesian;
    std::vector<mrcpp::GaussExp<3> *> orbitals;

    void readAOExpansion(Intgrl &intgrl);
    void transformToSpherical();
};

} //namespace gto_utils
} //namespace mrchem
