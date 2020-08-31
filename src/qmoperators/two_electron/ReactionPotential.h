#pragma once
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/two_electron/SCRF.h"

using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {
class ReactionPotential final : public QMPotential {
public:
    ReactionPotential(OrbitalVector_p Phi_p, SCRF help, bool var = false);
    ~ReactionPotential() override { free(NUMBER::Total); }
    friend class ReactionOperator;

    bool getRunVariational() const { return this->variational; }
    SCRF getHelper() { return this->helper; }
    double getNuclearEnergy() { return this->helper.getNuclearEnergy(); }
    double getElectronicEnergy() { return this->helper.getElectronicEnergy(); }
    double getTotalEnergy() { return this->helper.getTotalEnergy(); }
    void updateMOResidual(double const err_t) { this->helper.mo_residual = err_t; }

    QMFunction &getCurrentReactionPotential() { return this->helper.getCurrentReactionPotential(); }
    QMFunction &getPreviousReactionPotential() { return this->helper.getPreviousReactionPotential(); }
    QMFunction &getCurrentDifferenceReactionPotential() { return this->helper.getCurrentDifferenceReactionPotential(); }

    QMFunction &getCurrentGamma() { return this->helper.getCurrentGamma(); }
    QMFunction &getPreviousGamma() { return this->helper.getPreviousGamma(); }
    QMFunction &getCurrentDifferenceGamma() { return this->helper.getCurrentDifferenceGamma(); }

protected:
    void clear();

private:
    bool first_iteration = true;
    bool variational; // determines if the Reaction potential will be optimized in its own loop each SCF iteration or if
                      // it will converge together with the SCF procedure
    OrbitalVector_p Phi;
    SCRF helper;

    void setup(double prec);
};

} // namespace mrchem
