#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <vector>

#include "MRCPP/Gaussians"

class Intgrl;

class OrbitalExp {
public:
    OrbitalExp(Intgrl &intgrl);
    virtual ~OrbitalExp();

    int size() const { return orbitals.size(); }
    int getAngularMomentum(int n) const;

    mrcpp::GaussExp<3> &getOrbital(int n) { return *orbitals[n]; }

    void rotate(Eigen::MatrixXd &U);

protected:
    bool cartesian;
    std::vector<mrcpp::GaussExp<3> *> orbitals;

    void readAOExpansion(Intgrl &intgrl);
    void transformToSpherical();
};

