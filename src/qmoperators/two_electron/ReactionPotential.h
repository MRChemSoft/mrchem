#pragma once
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "chemistry/chemistry_fwd.h"
#include "qmfunctions/Density.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "scf_solver/KAIN.h"

using namespace mrcpp;

namespace mrchem {

class ReactionPotential final : public QMPotential {
public:
    ReactionPotential(mrcpp::PoissonOperator *P,
                      mrcpp::DerivativeOperator<3> *D,
                      Cavity *C,
                      const Nuclei &nucs,
                      OrbitalVector *Phi,
                      int hist,
                      double eps_i = 1.0,
                      double eps_o = 2.0,
                      bool islin = false);
    ~ReactionPotential() = default;

    double &getTotalEnergy();
    double &getElectronicEnergy();
    double &getNuclearEnergy();
    double &getElectronIn();
    QMFunction &getGamma() { return gamma; }
    QMFunction &getGammanp1() { return gammanp1; }
    void setGamma(QMFunction new_gamma) { this->gamma = new_gamma; }
    void setGammanp1(QMFunction new_gamma) { this->gammanp1 = new_gamma; }
    friend class ReactionOperator;

protected:
    void clear();

private:
    Cavity *cavity;
    Nuclei nuclei;
    OrbitalVector *orbitals;
    mrcpp::PoissonOperator *poisson;
    mrcpp::DerivativeOperator<3> *derivative;

    Density rho_tot;
    Density rho_el;
    Density rho_nuc;
    QMFunction cavity_func;

    QMFunction gamma;
    QMFunction gammanp1;

    int history;

    double d_coefficient = std::log(e_i / e_o);
    double electronicEnergy;
    double nuclearEnergy;
    double totalEnergy;
    double electronsIn;
    double e_i;
    double e_o;
    bool is_lin;

    void setRhoEff(QMFunction &rho_eff_func, std::function<double(const mrcpp::Coord<3> &r)> eps);
    void setGamma(QMFunction const &inv_eps_func,
                  QMFunction &gamma_func,
                  QMFunction &V_func,
                  mrcpp::FunctionTreeVector<3> &d_cavity);
    void accelerateConvergence(QMFunction &diff_func, QMFunction &temp, KAIN &kain);
    void setup(double prec);
};

} // namespace mrchem
