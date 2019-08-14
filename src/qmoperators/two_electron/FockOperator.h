#pragma once

#include "qmoperators/RankZeroTensorOperator.h"

/** @class FockOperator
 *
 * @brief Operator containing the standard SCF operators
 *
 * This is a simple collection of operators used in ground state SCF calculations.
 * The operator is separated into kinetic and potential parts, since the MW way of
 * solving the SCF equations is to invert the kinetic part, and apply the potential
 * part as usual.
 */

namespace mrchem {

class SCFEnergy;
class KineticOperator;
class NuclearOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;
class ElectricFieldOperator;
class ReactionOperator;

class FockOperator final : public RankZeroTensorOperator {
public:
    FockOperator(std::shared_ptr<KineticOperator> t = nullptr,
                 std::shared_ptr<NuclearOperator> v = nullptr,
                 std::shared_ptr<CoulombOperator> j = nullptr,
                 std::shared_ptr<ExchangeOperator> k = nullptr,
                 std::shared_ptr<XCOperator> xc = nullptr,
                 std::shared_ptr<ElectricFieldOperator> ext = nullptr,
                 std::shared_ptr<ReactionOperator> reo = nullptr);

    RankZeroTensorOperator &kinetic() { return this->T; }
    RankZeroTensorOperator &potential() { return this->V; }
    RankZeroTensorOperator &perturbation() { return this->H_1; }

    std::shared_ptr<KineticOperator> &getKineticOperator() { return this->kin; }
    std::shared_ptr<NuclearOperator> &getNuclearOperator() { return this->nuc; }
    std::shared_ptr<CoulombOperator> &getCoulombOperator() { return this->coul; }
    std::shared_ptr<ExchangeOperator> &getExchangeOperator() { return this->ex; }
    std::shared_ptr<XCOperator> &getXCOperator() { return this->xc; }
    std::shared_ptr<ElectricFieldOperator> &getExtOperator() { return this->ext; }
    std::shared_ptr<ReactionOperator> &getReactionOperator() { return this->Ro; }

    void rotate(const ComplexMatrix &U);

    void build(double exx = 1.0);
    void setup(double prec);
    void clear();

    SCFEnergy trace(OrbitalVector &Phi, const Nuclei &nucs);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

private:
    double exact_exchange{1.0};
    RankZeroTensorOperator T;   ///< Total kinetic energy operator
    RankZeroTensorOperator V;   ///< Total potential energy operator
    RankZeroTensorOperator H_1; ///< Perturbation operators

    std::shared_ptr<KineticOperator> kin;
    std::shared_ptr<NuclearOperator> nuc;
    std::shared_ptr<CoulombOperator> coul;
    std::shared_ptr<ExchangeOperator> ex;
    std::shared_ptr<XCOperator> xc;
    std::shared_ptr<ElectricFieldOperator> ext; ///< Total external potential
    std::shared_ptr<ReactionOperator> Ro;       ///< Reaction field operator
};

} // namespace mrchem
