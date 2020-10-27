#pragma once
#include "chemistry/Cavity.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmoperators/one_electron/QMPotential.h"
#include "qmoperators/two_electron/SCRF.h"

using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {
/** @brief class containing the solvent-substrate interaction reaction potential
 *
 */
class ReactionPotential final : public QMPotential {
public:
    ReactionPotential(OrbitalVector_p Phi_p, SCRF helper);
    ~ReactionPotential() override { free(NUMBER::Total); }
    friend class ReactionOperator;

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
    void setTesting() { this->first_iteration = false; }

protected:
    void clear();

private:
    bool first_iteration = true;
    OrbitalVector_p Phi;
    SCRF helper;

    void setup(double prec);
};

} // namespace mrchem
