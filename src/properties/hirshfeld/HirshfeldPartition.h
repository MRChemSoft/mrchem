#pragma once

#include <Eigen/Dense>
#include <mrchem.h>
#include "chemistry/Molecule.h"
#include <string>
#include "chemistry/Nucleus.h"
#include "properties/hirshfeld/HirshfeldInterpolator.h"

/**
 * @brief Class for computing the Hirshfeld partitioning of a molecule. Reads and interpolates the Hirshfeld partitioning data.
 * can create MW representation of the Hirshfeld partitioning functions.
 */
class HirshfeldPartition{

    public:
    /**
     * @brief Construct a new Hirshfeld Partition object
     * @param mol The molecule for which the Hirshfeld partitioning is to be computed
     * @param data_dir The directory containing the Hirshfeld partitioning data
     */
    HirshfeldPartition(const mrchem::Molecule &mol, std::string data_dir);

    /**
     * @brief Get the integral rho * w_i for a given atom i
     */
    double getHirshfeldPartitionIntegral(int index, mrcpp::ComplexFunction &rho, double prec) const;

    protected:

    /**
     * @brief Evaluate the analytic, interpolated Hirshfeld partitioning function at a given point
     */
    double evalf(const mrcpp::Coord<3> &r, int iAt) const;

    /**
     * @brief Evaluate the log sum of the exponential of the log atomic densities at a given point.
     */
    double lseLogDens(const mrcpp::Coord<3> &r) const;

    /**
     * @brief The nuclei of the molecule
     */
    std::shared_ptr<mrchem::Nuclei> nucs{nullptr};

    /**
     * @brief The number of nuclei in the molecule
     */
    int nNucs;

    /**
     * @brief The atomic density interpolators for the nuclei
     */
    std::vector<HirshfeldRadInterpolater> logDensities;


};