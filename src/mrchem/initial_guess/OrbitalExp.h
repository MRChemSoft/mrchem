#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <vector>

class Intgrl;
template<int D> class GaussExp;

class OrbitalExp {
public:
    OrbitalExp(Intgrl &intgrl);
    virtual ~OrbitalExp();

    int size() const { return orbitals.size(); }
    int getAngularMomentum(int n) const;

    GaussExp<3> &getOrbital(int n) { return *orbitals[n]; }

    void rotate(Eigen::MatrixXd &U);

protected:
    bool cartesian;
    std::vector<GaussExp<3> *> orbitals;

    void readAOExpansion(Intgrl &intgrl);
    void transformToSpherical();
};

