#pragma once
#include <mrchem.h>
#include "pseudopotential/sphericalHarmonics.h"

/**
 * @class ProjectorFunction
 * @brief Represents a nonlocal projector function for pseudopotential calculations.
 *
 * Each projector is characterized by its position, angular momentum quantum numbers (l, m),
 * a radial length parameter (rl), and a vector index (i). The projector is projected onto
 * the multiresolution basis upon construction.
 */
class ProjectorFunction{

public:
    /**
     * @brief Constructs a ProjectorFunction and projects it onto the multiresolution basis.
     * @param pos Position of the corresponding atom.
     * @param rl Radial length parameter of the projector.
     * @param i Index of the radial function p_i.
     * @param l Angular momentum quantum number.
     * @param m Magnetic quantum number.
     * @param prec Precision for the multiresolution projection.
     */
    ProjectorFunction(mrcpp::Coord<3> pos, double rl, int i, int l, int m, double prec);

    mrcpp::Coord<3> pos;                                ///< Position of the corresponding atom.
    int i;                                               ///< Index of the radial function p_i.
    int l;                                               ///< Angular momentum quantum number.
    int m;                                               ///< Magnetic quantum number.
    double rl;                                           ///< Radial length parameter.
    double prec;                                         ///< Precision used for the projection.
    std::shared_ptr<mrcpp::CompFunction<3>> projector_ptr; ///< Projected multiresolution representation.

private:
    /**
     * @brief Function pointer to the spherical harmonic for the given (l, m).
     */
    double (*s)(const std::array<double, 3> &r, const double &normr);

    /**
     * @brief Sets the spherical harmonic function pointer based on (l, m).
     * @param l Angular momentum quantum number.
     * @param m Magnetic quantum number.
     */
    void switch_sperics(int l, int m);

};
