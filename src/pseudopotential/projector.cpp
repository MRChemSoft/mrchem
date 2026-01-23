
#include "pseudopotential/projector.h"
#include <math.h>
#include <fstream>
#include <MRCPP/Printer>

#include <string>

// #include "mrchem.h"

/**
 * @brief Constructs a Projector object.
 *
 * This constructor initializes a Projector object with the given parameters.
 *
 * @param pos The position vector of the corresponding atom.
 * @param rl The radial length of the projector.
 * @param i Vector index of the projector.
 * @param l The angular momentum quantum number of the projector.
 * @param m The magnetic quantum number of the projector.
 * @param prec The precision of the projector.
 */
ProjectorFunction::ProjectorFunction(mrcpp::Coord<3> pos, double rl, int i, int l, int m, double prec) {
    this->pos = pos;
    this->rl = rl;
    this->i = i;
    int ii = i + 1;
    this->l = l;
    this->m = m;
    this->prec = 0.1 * prec;
    // select the spherical harmonic function based on the angular momentum and magnetic quantum number
    switch_sperics(l, m);
    double prefactor = std::sqrt(2.0) / (std::pow(rl, l + (4.0 * ii - 1) / 2.0) * std::sqrt(tgamma( l + (4.0 * ii - 1.0) / 2.0 )) );

    auto project_analytic = [this, prefactor, ii](const std::array<double, 3> &r) -> double {
        std::array<double, 3> rprime = {r[0] - this->pos[0], r[1] - this->pos[1], r[2] - this->pos[2]};
        double normr = std::sqrt( rprime[0] * rprime[0] + rprime[1] * rprime[1] + rprime[2] * rprime[2]);
        return prefactor * std::pow(normr, 2 * (ii - 1)) * std::exp(- 0.5 * normr * normr / (this->rl * this-> rl) ) * this->s(rprime, normr);
    };
    // auto op = (*this);
    // mrcpp::ComplexFunction f;
    projector_ptr = std::make_shared<mrcpp::CompFunction<3>>();

    double sigma = 0.6;
    auto gauss = [this, sigma](const std::array<double, 3> &r) -> double {
        std::array<double, 3> rprime = {r[0] - this->pos[0], r[1] - this->pos[1], r[2] - this->pos[2]};
        double normr = std::sqrt( rprime[0] * rprime[0] + rprime[1] * rprime[1] + rprime[2] * rprime[2]);
        double gaussNormalization = 1.0 / std::pow(2.0 * M_PI * sigma * sigma, 1.5);
        return std::exp(- 0.5 * normr * normr / (sigma * sigma) );
    };

    mrcpp::project(*projector_ptr, static_cast<std::function<double(const mrcpp::Coord<3>&)>>(gauss), prec);
    mrcpp::project(*projector_ptr,  static_cast<std::function<double(const mrcpp::Coord<3>&)>>(project_analytic), prec);
    // mrcpp::cplxfunc::deep_copy(op, f);

    double nrm = projector_ptr->norm();


    if (std::abs(nrm - 1.0) > 10 * prec) {
        MSG_ABORT("Projection of projector failed: norm=" << nrm << ", nodes=" << projector_ptr->getNNodes() << ", l=" << l << ", m=" << m << ", rl=" << rl << ", i=" << i);
    }

}

void ProjectorFunction::switch_sperics(int l, int m){
    this->s = get_spherical_harmonics(l, m);
}
