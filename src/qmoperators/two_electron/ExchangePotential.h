#pragma once

#include "QMOperator.h"
#include "qmfunctions.h"

namespace mrchem {

/** @class ExchangePotential
 *
 *  @brief Hartree-Fock exchange potential defined by a particular set of orbitals
 *
 * The operator is defined as the Hartree-Fock exchange arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. The internal exchange potentials (the operator applied to it's
 * own orbitals) can be precomputed and stored for fast retrieval. Option to use
 * screening based on previous calculations of the internal exchange (make sure that
 * the internal orbitals haven't been significantly changed since the last time the
 * operator was set up, e.g. through an orbital rotation).
 */

class ExchangePotential final : public QMOperator {
public:
    ExchangePotential(mrcpp::PoissonOperator &P, OrbitalVector &phi, bool s);
    ~ExchangePotential() { }

    void rotate(const ComplexMatrix &U);

    void setupInternal(double prec);

protected:
    bool screen;                ///< Apply screening in exchange evaluation
    DoubleVector tot_norms;     ///< Total norms for use in screening
    DoubleMatrix part_norms;    ///< Partial norms for use in screening
    OrbitalVector exchange;     ///< Precomputed exchange orbitals from the occupied orbital set

    // Pointers to external objects, ownership outside this class
    OrbitalVector *orbitals;         ///< Orbitals defining the exchange operator
    mrcpp::PoissonOperator *poisson; ///< Poisson operator to compute orbital contributions

    void setup(double prec);
    void clear();

    Orbital apply(Orbital phi_p);
    Orbital dagger(Orbital phi_p);

    using QMOperator::apply;
    using QMOperator::dagger;

    int testPreComputed(Orbital phi_p) const;
    double getScaledPrecision(int i, int j) const;
    double getSpinFactor(Orbital phi_i, Orbital phi_j) const;

    Orbital calcExchange(Orbital phi_p);

    void calcInternal(int i);
    void calcInternal(int i, int j);
    void calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j);
    void calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j, Orbital *V_ij );
};

} //namespace mrchem
