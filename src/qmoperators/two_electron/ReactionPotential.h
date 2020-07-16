#pragma once
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/two_electron/SCRF.h"

namespace mrchem {
class SCRF;
class ReactionPotential final : public QMPotential {
public:
    ReactionPotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                      std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                      int hist);
    ~ReactionPotential() override { free(NUMBER::Total); }
    friend class ReactionOperator;
    friend class SCRF;

    bool getRunVariational() const { return variational; }
    void setRunVariational(bool var) { this->variational = var; }
    std::shared_ptr<SCRF> getHelper() { return this->helper; }
    void updateTotalDensity(OrbitalVector Phi, double prec) { this->helper->updateTotalDensity(Phi, prec); }
    void setHelper(std::shared_ptr<SCRF> new_helper) { this->helper = new_helper; }
    double getNuclearEnergy() { return this->helper->getNuclearEnergy(); }
    double getElectronicEnergy() { return this->helper->getElectronicEnergy(); }
    double getTotalEnergy() { return this->helper->getTotalEnergy(); }
    void initial_setup() { this->run_once == true; }
    void updateMOResidual(double const err_t) { this->mo_residual = err_t; }

protected:
    void clear();

private:
    bool variational; // determines if the Reaction potential will be optimized in its own loop each SCF iteration or if
                      // it will converge together with the SCF procedure
    bool run_once;
    bool run_hybrid = false;
    bool run_absolute = false;
    int history;
    double mo_residual;
    std::shared_ptr<mrcpp::PoissonOperator> poisson;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;
    std::shared_ptr<SCRF> helper;

    void accelerateConvergence(QMFunction &diff_func, QMFunction &temp, KAIN &kain);
    void setup(double prec);
};

} // namespace mrchem
