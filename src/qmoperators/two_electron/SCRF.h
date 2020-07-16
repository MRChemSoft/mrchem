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
// using ReactionPotential_p = std::shared_ptr<mrchem::ReactionPotential>;

namespace mrchem {
class ReactionPotential;
class SCRF final {
public:
    SCRF(std::shared_ptr<ReactionPotential> Rp, Nuclei N, Permittivity e, OrbitalVector_p phi, double orb_prec);
    friend class ReactionPotential;
    void updateTotalDensity(OrbitalVector Phi,
                            double prec); // pass the electron orbitals and computes the total density
    void UpdateExternalDensity(Density new_density) { this->rho_ext = new_density; }
    void updateDifferencePotential(QMFunction diff_potential);
    void updateReactionPotential(QMFunction reac_potential);
    std::shared_ptr<mrchem::ReactionPotential> getReactionPotential() { return this->reaction_potential; }
    QMFunction &getDifferencePotential() { return this->difference_potential; }
    KAIN &getKain() { return this->reaction_optimizer; }

protected:
    void clear();

private:
    double apply_prec;
    Permittivity epsilon;
    Density rho_nuc;
    Density rho_ext;
    Density rho_tot;
    QMFunction difference_potential;
    std::shared_ptr<mrchem::ReactionPotential> reaction_potential;
    KAIN reaction_optimizer;
    mrcpp::FunctionTreeVector<3> d_cavity; // Vector containing the 3 partial derivatives of the cavity function

    QMFunctionVector makeTerms(DerivativeOperator_p derivative, PoissonOperator_p poisson, double prec);
    void updateGamma(QMFunction &gamma_func,
                     const DerivativeOperator_p derivative,
                     QMFunction potential_nm1,
                     const double prec);
    double getNuclearEnergy();
    double getElectronicEnergy();
    double getTotalEnergy();
    void resetQMFunction(QMFunction &function);
};
} // namespace mrchem
