#pragma once

#include "QMOperator.h"

namespace mrchem {

/** @class ExchangePotential
 *  @brief Hartree-Fock exchange potential for a given orbital.
 *  @author Stig Rune Jensen
 *  @date 2015, revised March 2018
 *
 */

class ExchangePotential final : public QMOperator {
public:
    ExchangePotential(mrcpp::PoissonOperator &P, OrbitalVector &phi, bool s);
    ~ExchangePotential() { }

    void rotate(const ComplexMatrix &U);

    void setup(double prec);
    void clear();

protected:
    bool screen;                ///< Apply screening in exchange evaluation
    DoubleVector tot_norms;     ///< Total norms for use in screening
    DoubleMatrix part_norms;    ///< Partial norms for use in screening
    OrbitalVector exchange;     ///< Precomputed exchange orbitals from the occupied orbital set

    // Pointers to external objects, ownership outside this class
    OrbitalVector *orbitals;         ///< Orbitals defining the exchange operator
    mrcpp::PoissonOperator *poisson; ///< Poisson operator to compute orbital contributions

    Orbital apply(Orbital phi_p);
    Orbital dagger(Orbital phi_p);

    using QMOperator::apply;
    using QMOperator::dagger;

    int testPreComputed(Orbital phi_p);
    double getScaledPrecision(int i, int j) const;
    double getSpinFactor(Orbital phi_i, Orbital phi_j) const;

    Orbital calcExchange(Orbital phi_p);

    void calcInternalExchange();
    void calcInternal(int i);
    void calcInternal(int i, int j);
    void calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j);
    void calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j, Orbital *V_ij );
};

} //namespace mrchem
