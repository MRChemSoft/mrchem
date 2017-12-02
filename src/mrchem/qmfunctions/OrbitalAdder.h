#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "MWAdder.h"
#include "GridGenerator.h"
#include "constants.h"

class Orbital;
class OrbitalVector;

class OrbitalAdder {
public:
    OrbitalAdder(double prec, int max_scale, int work_vec_max = workOrbVecSize);
    virtual ~OrbitalAdder() { }

    void setPrecision(double prec) { this->add.setPrecision(prec); }

    void operator()(Orbital &phi_ab,
                    std::complex<double> a, Orbital &phi_a,
                    std::complex<double> b, Orbital &phi_b,
                    bool union_grid);

    void operator()(Orbital &out,
                    std::vector<std::complex<double> > &coefs,
                    std::vector<Orbital *> &inp,
                    bool union_grid);

    void operator()(OrbitalVector &out,
                    double a, OrbitalVector &inp_a,
                    double b, OrbitalVector &inp_b,
                    bool union_grid);

    void operator()(Orbital &out,
                    const Eigen::VectorXd &c,
                    OrbitalVector &inp,
                    bool union_grid);

    void rotate(OrbitalVector &out, const Eigen::MatrixXd &U, OrbitalVector &inp);
    void rotate_P(OrbitalVector &out, const Eigen::MatrixXd &U, OrbitalVector &inp);
    void rotate(OrbitalVector &out, const Eigen::MatrixXd &U);

    void inPlace(Orbital &out, std::complex<double> c, Orbital &inp);
    void inPlace(Orbital &out, std::vector<std::complex<double> > &coefs, std::vector<Orbital *> &inp, bool union_grid);
    void inPlace(Orbital &out, const Eigen::VectorXd &c, OrbitalVector &inp, bool union_grid);
    void inPlace(OrbitalVector &out, double c, OrbitalVector &inp);

    void orthogonalize(OrbitalVector &out);
    void orthogonalize(Orbital &out, Orbital &inp);
    void orthogonalize(Orbital &out, OrbitalVector &inp);
    void orthogonalize(OrbitalVector &out, Orbital &inp);
    void orthogonalize(OrbitalVector &out, OrbitalVector &inp);

protected:
    int workVecMax; //max number of orbitals to store temporarily while sending
    MWAdder<3> add;
    GridGenerator<3> grid;
};

