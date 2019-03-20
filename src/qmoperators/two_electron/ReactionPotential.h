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

    int history;

    double d_coefficient = std::log(e_i / e_o);
    double electronicEnergy;
    double nuclearEnergy;
    double totalEnergy;
    double e_i;
    double e_o;
    bool is_lin;

    void setRhoEff(QMFunction &rho_eff_func, std::function<double(const mrcpp::Coord<3> &r)> eps);
    void setGamma(QMFunction const &inv_eps_func,
                  QMFunction &gamma_func,
                  QMFunction &temp_func1,
                  mrcpp::FunctionTreeVector<3> &d_cavity);
    void grad_G(QMFunction &gamma_func,
                QMFunction &rho_tot,
                QMFunction &grad_G_func,
                std::function<double(const mrcpp::Coord<3> &r)> eps);
    void accelerateConvergence(QMFunction &diff_func, QMFunction &temp, KAIN &kain);
    void setup(double prec);
};

} // namespace mrchem
