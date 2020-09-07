#pragma once

#include "chemistry/Nucleus.h"
#include "chemistry/Permittivity.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/QMFunction.h"
#include "qmfunctions/qmfunction_fwd.h"
#include "scf_solver/KAIN.h"

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {
class ReactionPotential;
class SCRF final {
public:
    SCRF(Permittivity e,
         const Nuclei &N,
         PoissonOperator_p P,
         DerivativeOperator_p D,
         double orb_prec,
         int kain_hist = 0,
         int max_iter = 100,
         bool accelerate_Vr = true,
         bool run_hybrid = true,
         bool run_absolute = false);
    friend class ReactionPotential;
    void UpdateExternalDensity(Density new_density) { this->rho_ext = new_density; }

    QMFunction &getCurrentReactionPotential() { return this->Vr_n; }
    QMFunction &getPreviousReactionPotential() { return this->Vr_nm1; }
    QMFunction &getCurrentDifferenceReactionPotential() { return this->dVr_n; }

    QMFunction &getCurrentGamma() { return this->gamma_n; }
    QMFunction &getPreviousGamma() { return this->gamma_nm1; }
    QMFunction &getCurrentDifferenceGamma() { return this->dgamma_n; }

    void updateMOResidual(double const err_t) { this->mo_residual = err_t; }

protected:
    void clear();

private:
    bool accelerate_Vr;
    bool run_hybrid;
    bool run_absolute;

    int max_iter;
    int history;
    double apply_prec;
    double mo_residual;

    Permittivity epsilon;

    Density rho_nuc;
    Density rho_ext;
    Density rho_tot;

    QMFunction Vr_n;
    QMFunction dVr_n;
    QMFunction Vr_nm1;

    QMFunction gamma_n;
    QMFunction dgamma_n;
    QMFunction gamma_nm1;

    mrcpp::FunctionTreeVector<3> d_cavity; // Vector containing the 3 partial derivatives of the cavity function
    DerivativeOperator_p derivative;
    PoissonOperator_p poisson;

    void setDCavity();

    void computeDensities(OrbitalVector &Phi);
    void computeGamma(QMFunction Potential, QMFunction &out_gamma);

    QMFunction solvePoissonEquation(const QMFunction &ingamma);

    void accelerateConvergence(QMFunction &dfunc, QMFunction &func, KAIN &kain);

    void variationalSCRF(QMFunction V_vac);
    void nestedSCRF(QMFunction V_vac);
    QMFunction &setup(double prec, const OrbitalVector_p &Phi, bool variational = false);

    double getNuclearEnergy();
    double getElectronicEnergy();
    double getTotalEnergy();
    void resetQMFunction(QMFunction &function);
    void updateCurrentReactionPotential(QMFunction &Vr_np1);
    void updateCurrentGamma(QMFunction &gamma_np1);
};
} // namespace mrchem
