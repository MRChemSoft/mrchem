#pragma once
#include "chemistry/Nucleus.h"
#include <vector>
#include "utils/PolyInterpolator.h"
#include "qmoperators/QMPotential.h"


namespace mrchem {

/** @class AZora
 *
 * @brief Operator defining the AZora potential based on molecular data.
 *
 * Inherits from QMPotential and adds functionality to utilize an mrchem::Molecule
 * for constructing the potential.
 *
 */
class AZoraPotential : public QMPotential {
public:
    /**
     * Constructor that takes a molecule and initializes the azora potential.
     * @param molecule The molecule used to construct the potential.
     * @param adap Adaptive parameter from QMPotential.
     * @param shared Determines if the base potential is shared.
     */
    AZoraPotential(Nuclei nucs, int adap, std::string azora_dir, bool shared = false, double c = 137.035999084) 
        : QMPotential(adap, shared) {
        this->nucs = nucs;
        this->c = c;
        this->azora_dir = azora_dir;
        initAzoraPotential();
    }

    /**
     * Copy constructor.
     * @param other The other instance to copy from.
     */
    AZoraPotential(const AZoraPotential& other) 
        : QMPotential(other) {
        this->nucs = other.nucs;
        this->prec = other.prec;
        this->c = other.c;
        this->azora_dir = other.azora_dir;
        initAzoraPotential();
    }

    /**
     * Destructor.
     */
    virtual ~AZoraPotential(){
        free(mrchem::NUMBER::Total);
        isProjected = false;
    }

    // Delete copy assignment to prevent copying
    AZoraPotential& operator=(const AZoraPotential&) = delete;

    /**
     * Project the potential from the analytic function stored in this object.
     * @param prec Projection precision.
     */
    void project(double prec);

protected:
    Nuclei nucs; // The nuclei of the molecule
    double prec; // The precision parameter
    double c;    // The speed of light
    std::string azora_dir; // The directory containing the azora potential data
    std::vector<interpolation_utils::PolyInterpolator> atomicPotentials;
    bool isProjected = false;

    double evalf_analytic(const mrcpp::Coord<3>& r);

    /**
     * Initialize the azora potential based on the molecule.
     * This method would typically setup the real and imaginary function trees
     * representing the potential.
     */
    void initAzoraPotential();
};

} // namespace mrchem