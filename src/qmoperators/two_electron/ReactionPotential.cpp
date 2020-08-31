#include "ReactionPotential.h"
#include "MRCPP/MWOperators"
#include "MRCPP/Plotter"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "scf_solver/KAIN.h"
#include "utils/print_utils.h"
#include <string>

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using Cavity_p = std::shared_ptr<mrchem::Cavity>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

ReactionPotential::ReactionPotential(OrbitalVector_p Phi_p, SCRF help, bool var)
        : QMPotential(1, false)
        , variational(var)
        , Phi(Phi_p)
        , helper(help) {}

void ReactionPotential::setup(double prec) {
    setApplyPrec(prec);
    QMFunction &temp = *this;
    if (this->first_iteration) {
        this->first_iteration = false;
        return;
    }
    qmfunction::deep_copy(temp, this->helper.setup(prec, this->Phi, this->variational));
}

void ReactionPotential::clear() {
    //    QMFunction::free(NUMBER::Real);
    clearApplyPrec();
    this->helper.clear();
}

} // namespace mrchem
