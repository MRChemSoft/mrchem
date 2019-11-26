#pragma once
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "chemistry/chemistry_fwd.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "scf_solver/KAIN.h"

namespace mrchem {

class ReactionPotential final : public QMPotential {
public:
    ReactionPotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                      std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                      std::shared_ptr<mrchem::Cavity> C,
                      const Nuclei &nucs,
                      std::shared_ptr<mrchem::OrbitalVector> Phi,
                      int hist,
                      double eps_i,
                      double eps_o,
                      bool islin);
    ~ReactionPotential() = default;

    friend class ReactionOperator;
    double &getTotalEnergy();
    double &getElectronicEnergy();
    double &getNuclearEnergy();
    double &getElectronIn();
    QMFunction &getGamma() { return gamma; }
    QMFunction &getGammanp1() { return gammanp1; }
    bool &getRunVariational() { return variational; }
    void setGamma(QMFunction new_gamma) { this->gamma = new_gamma; }
    void setGammanp1(QMFunction new_gamma) { this->gammanp1 = new_gamma; }
    void setRunVariational(bool var) { this->variational = var; }

protected:
    void clear();

private:
    std::shared_ptr<Cavity> cavity;
    Nuclei nuclei;
    std::shared_ptr<OrbitalVector> orbitals;
    std::shared_ptr<mrcpp::PoissonOperator> poisson;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;

    Density rho_tot;
    Density rho_el;
    Density rho_nuc;
    QMFunction cavity_func;

    mrcpp::FunctionTreeVector<3> d_cavity;

    QMFunction gamma;
    QMFunction gammanp1;

    int history;

    double d_coefficient; // factor with which rescale the derivative of the cavity function
    double electronicEnergy{0.0};
    double nuclearEnergy{0.0};
    double totalEnergy{0.0};
    double electronsIn{0.0};
    double e_i;       // permitivity of free space, 1 inside the cavity
    double e_o;       // dielectric constant characteristic of the solvent
    bool is_lin;      // determines if the dielectric function will be implemented linearly or exponentially
    bool variational; // determines if the Reaction potential will be optimized in its own loop each SCF iteration or if
                      // it will converge together with the SCF procedure

    void setRhoEff(QMFunction &rho_eff_func, std::function<double(const mrcpp::Coord<3> &r)> eps);
    void setGamma(QMFunction const &inv_eps_func, QMFunction &gamma_func, QMFunction &V_func);
    void accelerateConvergence(QMFunction &diff_func, QMFunction &temp, KAIN &kain);
    void poissonSolver(QMFunction rho_eff_func, QMFunction *diff_func, QMFunction *V_np1_func, double *error);
    void SCRF(QMFunction *V_tot_func,
              QMFunction *V_vac_func,
              QMFunction *rho_eff_func,
              QMFunction &temp,
              const QMFunction &inv_eps_func);
    void setup(double prec);
};

} // namespace mrchem
